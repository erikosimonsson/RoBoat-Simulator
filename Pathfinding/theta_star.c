#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

#define INF DBL_MAX

typedef struct {
    int x, y;
} Point;

typedef struct {
    int x, y;
    double g, h, f;
    int px, py;
    bool opened, closed;
    int heap_index;
} Node;

typedef struct {
    Node **data;
    int size;
    int capacity;
} MinHeap;

static void heap_swap(MinHeap *heap, int i, int j);
static void heap_sift_up(MinHeap *heap, int i);
static void heap_sift_down(MinHeap *heap, int i);
static bool heap_init(MinHeap *heap, int capacity);
static void heap_free(MinHeap *heap);
static void heap_push(MinHeap *heap, Node *node);
static Node *heap_pop_min(MinHeap *heap);
static void heap_decrease_key(MinHeap *heap, Node *node);

static inline int grid_at(const int *grid, int width, int x, int y) {
    return grid[y * width + x];
}

static inline Node *node_at(Node *nodes, int width, int x, int y) {
    return &nodes[y * width + x];
}

static bool in_bounds(int x, int y, int width, int height) {
    return x >= 0 && x < width && y >= 0 && y < height;
}

static bool is_blocked(const int *grid, int width, int height, int x, int y) {
    return !in_bounds(x, y, width, height) || grid_at(grid, width, x, y) != 0;
}

static double distance_between(int x1, int y1, int x2, int y2) {
    double dx = (double)(x2 - x1);
    double dy = (double)(y2 - y1);
    return sqrt(dx * dx + dy * dy);
}

static double heuristic(int x1, int y1, int x2, int y2) {
    return distance_between(x1, y1, x2, y2);
}

static bool can_step(const int *grid, int width, int height, int x1, int y1, int x2, int y2) {
    if (!in_bounds(x2, y2, width, height) || is_blocked(grid, width, height, x2, y2)) {
        return false;
    }
    int dx = x2 - x1;
    int dy = y2 - y1;
    if (dx != 0 && dy != 0) {
        if (is_blocked(grid, width, height, x1 + dx, y1) || is_blocked(grid, width, height, x1, y1 + dy)) {
            return false;
        }
    }
    return true;
}

static bool line_of_sight(const int *grid, int width, int height, int x1, int y1, int x2, int y2) {
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x2 > x1) ? 1 : -1;
    int sy = (y2 > y1) ? 1 : -1;
    int err = dx - dy;
    int x = x1, y = y1;

    while (x != x2 || y != y2) {
        int prev_x = x, prev_y = y;
        int err2 = err * 2;
        if (err2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (err2 < dx) {
            err += dx;
            y += sy;
        }
        if (!can_step(grid, width, height, prev_x, prev_y, x, y)) {
            return false;
        }
    }
    return true;
}

static void init_nodes(Node *nodes, int width, int height, int goal_x, int goal_y) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Node *node = node_at(nodes, width, x, y);
            node->x = x;
            node->y = y;
            node->g = INF;
            node->h = heuristic(x, y, goal_x, goal_y);
            node->f = INF;
            node->px = x;
            node->py = y;
            node->opened = false;
            node->closed = false;
            node->heap_index = -1;
        }
    }
}

static void update_vertex(const int *grid, int width, int height, Node *nodes, MinHeap *open_heap, Node *current, Node *neighbor) {
    Node *parent = node_at(nodes, width, current->px, current->py);

    if (line_of_sight(grid, width, height, parent->x, parent->y, neighbor->x, neighbor->y)) {
        double new_g = parent->g + distance_between(parent->x, parent->y, neighbor->x, neighbor->y);
        if (new_g < neighbor->g) {
            neighbor->g = new_g;
            neighbor->px = parent->x;
            neighbor->py = parent->y;
            neighbor->f = neighbor->g + neighbor->h;
            if (!neighbor->opened) {
                neighbor->opened = true;
                heap_push(open_heap, neighbor);
            }
            else {
                heap_decrease_key(open_heap, neighbor);
            }
        }
    }
    else {
        double new_g = current->g + distance_between(current->x, current->y, neighbor->x, neighbor->y);
        if (new_g < neighbor->g) {
            neighbor->g = new_g;
            neighbor->px = current->x;
            neighbor->py = current->y;
            neighbor->f = neighbor->g + neighbor->h;
            if (!neighbor->opened) {
                neighbor->opened = true;
                heap_push(open_heap, neighbor);
            }
            else {
                heap_decrease_key(open_heap, neighbor);
            }
        }
    }
}

static bool reconstruct_path(Node *nodes, int width, int height, int start_x, int start_y, int goal_x, int goal_y, Point *path, int *path_length) {
    int max_path = width * height;
    int count = 0;
    int cx = goal_x, cy = goal_y;
    while (!(cx == start_x && cy == start_y)) {
        if (count >= max_path) {
            return false;
        }
        path[count++] = (Point){cx, cy};
        Node *node = node_at(nodes, width, cx, cy);
        if (node->px == cx && node->py == cy) {
            return false;
        }
        cx = node->px;
        cy = node->py;
    }
    path[count++] = (Point){start_x, start_y};
    for (int i = 0; i < count / 2; i++) {
        Point tmp = path[i];
        path[i] = path[count - 1 - i];
        path[count - 1 - i] = tmp;
    }
    *path_length = count;
    (void)height;
    return true;
}

bool theta_star(const int *grid, int width, int height, int start_x, int start_y, int goal_x, int goal_y, Point *path, int *path_length, double *path_cost) {
    static const int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    static const int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

    if (!grid || !path || !path_length) {
        return false;
    }
    if (!in_bounds(start_x, start_y, width, height) || !in_bounds(goal_x, goal_y, width, height)) {
        return false;
    }
    if (is_blocked(grid, width, height, start_x, start_y) || is_blocked(grid, width, height, goal_x, goal_y)) {
        return false;
    }
    Node *nodes = malloc((size_t)width * height * sizeof(Node));
    if (!nodes) {
        return false;
    }
    MinHeap open_heap;
    if (!heap_init(&open_heap, width * height)) {
        free(nodes);
        return false;
    }
    init_nodes(nodes, width, height, goal_x, goal_y);
    Node *start_node = node_at(nodes, width, start_x, start_y);
    start_node->g = 0.0;
    start_node->f = start_node->h;
    start_node->px = start_x;
    start_node->py = start_y;
    start_node->opened = true;
    heap_push(&open_heap, start_node);

    while (true) {
        Node *current = heap_pop_min(&open_heap);
        if (current == NULL) {
            heap_free(&open_heap);
            free(nodes);
            return false;
        }
        if (current->closed) {
            continue;
        }
        if (current->x == goal_x && current->y == goal_y) {
            bool result = reconstruct_path(nodes, width, height, start_x, start_y, goal_x, goal_y, path, path_length);
            if (result && path_cost) {
                *path_cost = current->g;
            }
            heap_free(&open_heap);
            free(nodes);
            return result;
        }
        current->closed = true;

        for (int i = 0; i < 8; i++) {
            int nx = current->x + dx[i];
            int ny = current->y + dy[i];
            if (!can_step(grid, width, height, current->x, current->y, nx, ny)) {
                continue;
            }
            Node *neighbor = node_at(nodes, width, nx, ny);
            if (neighbor->closed) {
                continue;
            }
            update_vertex(grid, width, height, nodes, &open_heap, current, neighbor);
        }
    }
}

static void draw_segment(char *display, int width, int height, int x1, int y1, int x2, int y2) {
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;
    int x = x1, y = y1;

    while (true) {
        if (in_bounds(x, y, width, height) && display[y * width + x] == '.') {
            display[y * width + x] = '+';
        }
        if (x == x2 && y == y2) {
            break;
        }
        int err2 = err * 2;
        if (err2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (err2 < dx) {
            err += dx;
            y += sy;
        }
    }
}

void print_grid_with_path(const int *grid, int width, int height, const Point *path, int path_length, int start_x, int start_y, int goal_x, int goal_y) {
    FILE *out = fopen("out.txt", "w");
    if (!out) {
        fprintf(stderr, "Error: could not open output file\n");
        return;
    }
    char *display = malloc((size_t)width * height);
    if (!display) {
        fprintf(stderr, "Error: memory allocation failed\n");
        fclose(out);
        return;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            display[y * width + x] = (grid_at(grid, width, x, y) == 0) ? '.' : '#';
        }
    }

    for (int i = 0; i < path_length - 1; i++) {
        draw_segment(display, width, height, path[i].x, path[i].y, path[i + 1].x, path[i + 1].y);
    }

    for (int i = 0; i < path_length; i++) {
        int x = path[i].x;
        int y = path[i].y;
        if (in_bounds(x, y, width, height)) {
            display[y * width + x] = 'W';
        }
    }
    display[start_y * width + start_x] = 'S';
    display[goal_y * width + goal_x] = 'G';

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            fputc(display[y * width + x], out);
        }
        fputc('\n', out);
    }
    free(display);
    fclose(out);
}

static bool load_grid_from_file(const char *filename, int **grid, int *width, int *height, int *start_x, int *start_y, int *goal_x, int *goal_y) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: could not open file %s\n", filename);
        return false;
    }

    char line[4096];
    int w = -1;
    int h = 0;
    int capacity = 0;
    int *g = NULL;
    bool found_start = false;
    bool found_goal = false;

    while (fgets(line, sizeof(line), file)) {
        size_t len = strcspn(line, "\r\n");
        line[len] = '\0';
        if (len == 0) {
            continue;
        }
        if (w == -1) {
            w = (int)len;
            if (w <= 0) {
                fprintf(stderr, "Error: invalid map width\n");
                fclose(file);
                return false;
            }
        } 
        else if ((int)len != w) {
            fprintf(stderr, "Error: inconsistent row width at row %d\n", h);
            free(g);
            fclose(file);
            return false;
        }
        if ((h + 1) * w > capacity) {
            int new_capacity = (capacity == 0) ? (w * 8) : (capacity * 2);
            while (new_capacity < (h + 1) * w) {
                new_capacity *= 2;
            }
            int *new_grid = realloc(g, (size_t)new_capacity * sizeof(int));
            if (!new_grid) {
                fprintf(stderr, "Error: memory allocation failed\n");
                free(g);
                fclose(file);
                return false;
            }
            g = new_grid;
            capacity = new_capacity;
        }

        for (int x = 0; x < w; x++) {
            char c = line[x];
            switch (c) {
                case '.':
                    g[h * w + x] = 0;
                    break;
                case '#':
                    g[h * w + x] = 1;
                    break;
                case 'S':
                    if (found_start) {
                        fprintf(stderr, "Error: multiple S cells found\n");
                        free(g);
                        fclose(file);
                        return false;
                    }
                    found_start = true;
                    *start_x = x;
                    *start_y = h;
                    g[h * w + x] = 0;
                    break;
                case 'G':
                    if (found_goal) {
                        fprintf(stderr, "Error: multiple G cells found\n");
                        free(g);
                        fclose(file);
                        return false;
                    }
                    found_goal = true;
                    *goal_x = x;
                    *goal_y = h;
                    g[h * w + x] = 0;
                    break;
                default:
                    fprintf(stderr, "Error: invalid character '%c' at row %d, col %d\n", c, h, x);
                    free(g);
                    fclose(file);
                    return false;
            }
        }
        h++;
    }
    fclose(file);
    if (w <= 0 || h <= 0) {
        fprintf(stderr, "Error: empty map file\n");
        free(g);
        return false;
    }
    if (!found_start) {
        fprintf(stderr, "Error: no S found in map\n");
        free(g);
        return false;
    }
    if (!found_goal) {
        fprintf(stderr, "Error: no G found in map\n");
        free(g);
        return false;
    }
    *grid = g;
    *width = w;
    *height = h;
    return true;
}

int main(int argc, char **argv) {
    const char *filename = argv[1];
    int *grid = NULL;
    int width = 0, height = 0;
    int start_x = -1, start_y = -1;
    int goal_x = -1, goal_y = -1;
    if (!load_grid_from_file(
            filename,
            &grid,
            &width,
            &height,
            &start_x,
            &start_y,
            &goal_x,
            &goal_y)) {
        return 1;
    }
    Point *path = malloc((size_t)width * height * sizeof(Point));
    if (!path) {
        fprintf(stderr, "Error: memory allocation failed.\n");
        free(grid);
        return 1;
    }
    int path_length = 0;
    double path_cost = 0.0;
    bool path_found = theta_star(
        grid, width, height,
        start_x, start_y,
        goal_x, goal_y,
        path, &path_length, &path_cost
    );

    if (path_found) {
        print_grid_with_path(grid, width, height, path, path_length, start_x, start_y, goal_x, goal_y);
    } else {
        printf("No path found\n");
    }
    free(path);
    free(grid);
    return 0;
}

//#------------------- MinHeap helper functions ------------------#

static void heap_swap(MinHeap *heap, int i, int j) {
    Node *tmp = heap->data[i];
    heap->data[i] = heap->data[j];
    heap->data[j] = tmp;

    heap->data[i]->heap_index = i;
    heap->data[j]->heap_index = j;
}

static void heap_sift_up(MinHeap *heap, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (heap->data[p]->f <= heap->data[i]->f) {
            break;
        }
        heap_swap(heap, i, p);
        i = p;
    }
}

static void heap_sift_down(MinHeap *heap, int i) {
    while (true) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = i;

        if (left < heap->size && heap->data[left]->f < heap->data[smallest]->f) {
            smallest = left;
        }
        if (right < heap->size && heap->data[right]->f < heap->data[smallest]->f) {
            smallest = right;
        }
        if (smallest == i) {
            break;
        }

        heap_swap(heap, i, smallest);
        i = smallest;
    }
}

static bool heap_init(MinHeap *heap, int capacity) {
    heap->data = malloc((size_t)capacity * sizeof(Node *));
    if (!heap->data) {
        return false;
    }
    heap->size = 0;
    heap->capacity = capacity;
    return true;
}

static void heap_free(MinHeap *heap) {
    free(heap->data);
    heap->data = NULL;
    heap->size = 0;
    heap->capacity = 0;
}

static void heap_push(MinHeap *heap, Node *node) {
    if (heap->size >= heap->capacity) {
        return;
    }
    heap->data[heap->size] = node;
    node->heap_index = heap->size;
    heap->size++;
    heap_sift_up(heap, node->heap_index);
}

static Node *heap_pop_min(MinHeap *heap) {
    if (heap->size == 0) {
        return NULL;
    }

    Node *min = heap->data[0];
    heap->size--;

    if (heap->size > 0) {
        heap->data[0] = heap->data[heap->size];
        heap->data[0]->heap_index = 0;
        heap_sift_down(heap, 0);
    }

    min->heap_index = -1;
    return min;
}

static void heap_decrease_key(MinHeap *heap, Node *node) {
    if (node->heap_index >= 0) {
        heap_sift_up(heap, node->heap_index);
    }
}