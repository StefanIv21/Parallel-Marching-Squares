// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
//strcutura trimisa ca argument la functia de start a thread-urilor
struct thread_arg {
	int id; //id-ul thread-ului
	int P; //numarul de thread-uri
	int p; //numarul de linii din grid
	int q; //numarul de coloane din grid
	unsigned char **grid; //grid-ul
	ppm_image *image; //imaginea 
    int step_x; //distanta pe x drintre punctele din grid
    int step_y; //distanta pe y drintre punctele din grid
    unsigned char sigma; //valoarea de threshold pentru grid ul binar
    ppm_image **contour_map; //map-ul de contur
    ppm_image * old_image;  //imaginea initiala
    pthread_barrier_t *barrier; //bariera
};

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

//functie pentru calculare minim
int min (int a ,int b)
{
    if(a>b)return b;
    return a;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}


void *threads_start(void *arg) {
    struct thread_arg* args = (struct thread_arg*) arg;

    int p = args->p;
    int q = args->q;
    int step_x = args->step_x;
    int step_y = args->step_y;
    unsigned char **grid = args->grid;
	ppm_image *image = args->image;
    ppm_image *old_image = args->old_image;
    unsigned char sigma = args->sigma;
    ppm_image **contour_map = args->contour_map;

    //start si end pentru algoritmul marching squares (in functie de p)
    int start = args->id * (double)p / args->P;
    int end = min((args->id + 1) * (double)p / args->P, p);

    //verific daca trebuie sa fac rescale(interpolare bicubica)
    if (old_image->x > RESCALE_X && old_image->y > RESCALE_Y) {

        //start si end pentru rescale
        int start1 = args->id * (double)image->x / args->P;
        int end1 = min((args->id + 1) * (double)image->x / args->P, image->x);

        uint8_t sample[3] ={0};
        for (int i = start1; i < end1; i++) {
            for (int j = 0; j < image->y; j++) {
            
                float u = (float)i / (float)(image->x - 1);
                float v = (float)j / (float)(image->y - 1);
                sample_bicubic(old_image, u, v, sample);

                image->data[i * image->y + j].red = sample[0];
                image->data[i * image->y + j].green = sample[1];
                image->data[i * image->y + j].blue = sample[2];
            }
        }
    }

    pthread_barrier_wait (args->barrier);

    // Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
    // Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
    // pixel values compare to the `sigma` reference value. The points are taken at equal distances
    // in the original image, based on the `step_x` and `step_y` arguments.

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;
 
    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them

    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    for (int j = 0; j < q; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    pthread_barrier_wait (args->barrier);

    // Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
    // type of contour which corresponds to each subgrid. It determines the binary value of each
    // sample fragment of the original image and replaces the pixels in the original image with
    // the pixels of the corresponding contour image accordingly.

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }

    return NULL;
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x,
                    struct thread_arg *arguments, pthread_t *threads, ppm_image *old_image) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);
    free(arguments);
    free(threads);

     if (old_image->x > RESCALE_X && old_image->y > RESCALE_Y) {
        free(old_image->data);
        free(old_image);
        free(image->data);
        free(image);
    }
    else
    {
        free(old_image->data);
        free(old_image);
    }

}

ppm_image *rescale_image(ppm_image *image) {
    

    // we only rescale downwards
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    // alloc memory for image
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    return new_image;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file>\n");
        return 1;
    }
    //citesc numarul de thread-uri
    int P = atoi(argv[3]);
    //status pentru terminarea thread-urilor
    void *status;
    int r;
    //alocare memorie pentru strcutura trimisa ca argument la functia de start a thread-urilor
    struct thread_arg *arguments;
    arguments = (struct thread_arg*) malloc(P * sizeof(struct thread_arg));
    if (!arguments) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    //alocare memorie pentru thread-uri
    pthread_t *threads;
    threads = (pthread_t*) malloc(P * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    //citire imagine
    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    //initializare bariera
    pthread_barrier_t   barrier;
    r = pthread_barrier_init (&barrier, NULL, P);
    if (r) {
        fprintf(stderr, "Unable to init barrier\n");
        exit(1);
    }

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();
   
    // 1. Rescale the image
    ppm_image *scaled_image = rescale_image(image);

    //// Builds a p x q grid of points
    int p = scaled_image->x / step_x;
    int q = scaled_image->y / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    //creare thread-uri si initializare argumente
    for (int i = 0; i < P; i++) {
		arguments[i].id = i;
		arguments[i].P = P;
		arguments[i].p =p;
        arguments[i].q =q;
        arguments[i].grid = grid;
        arguments[i].image=scaled_image;
        arguments[i].step_x= step_x;
        arguments[i].step_y =step_y;
        arguments[i].sigma=SIGMA;
        arguments[i].contour_map = contour_map;
        arguments[i].grid = grid;
        arguments[i].old_image = image;
        arguments[i].barrier = &barrier;

        r = pthread_create(&threads[i], NULL, threads_start, &arguments[i]);

		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}

    //asteptare thread-uri
	for (int i = 0; i < P; i++) {
		r = pthread_join(threads[i], &status);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}
    //distrugere bariera
    r = pthread_barrier_destroy (&barrier);
    if (r) {
        fprintf(stderr, "Unable to destroy barrier\n");
        exit(1);
    }

    // 4. Write output
    write_ppm(scaled_image, argv[2]);
    
    // 5. Free resources
    free_resources(scaled_image, contour_map, grid, step_x,arguments,threads,image);
    return 0;
}