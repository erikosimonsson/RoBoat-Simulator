#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#define MAX_ENVIRON_HEIGHT 10
#define MAX_ENVIRON_WIDTH 10
#define INF DBL_MAX

typedef struct Node {
    int x, y;
    double g, h, f;
    int px, py;
    bool closed, opened;
} Node;

typedef struct {
    int x, y;
} Point;

static int grid[MAX_ENVIRON_HEIGHT][MAX_ENVIRON_WIDTH] = {
    {0,0,0,0,0,0,0,0,0,0},
    {0,1,1,1,0,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,0,0},
    {0,0,0,1,0,0,1,1,1,0},
    {0,0,0,0,0,0,0,0,1,0},
    {0,0,0,0,1,1,1,0,1,0},
    {0,0,0,0,0,0,0,0,1,0},
    {0,0,1,1,1,1,0,0,0,0},
    {0,0,0,0,0,0,0,1,1,0},
    {0,0,0,0,0,0,0,0,0,0}
};

static Node nodes[MAX_ENVIRON_HEIGHT][MAX_ENVIRON_WIDTH];

static const int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static const int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

bool in_bounds(int x, int y) {
    return x >= 0 && x < MAX_ENVIRON_WIDTH && y >= 0 && y < MAX_ENVIRON_HEIGHT;
}

bool is_blocked(int x, int y) {
    return !in_bounds(x, y) || grid[y][x] != 0;
}

double heuristic(int x1, int y1, int x2, int y2) {
    double dx = (double)(x2 - x1);
    double dy = (double)(y2 - y1);
    return sqrt(dx * dx + dy * dy);
}

double distance_between(int x1, int y1, int x2, int y2) {
    return heuristic(x1, y1, x2, y2);
}

/**
 * Bresenham style traversal
 * Returns true if the segment from (x0, y0) to (x1, y1) does not cross blocked cells
 */
bool line_of_sight(int x0, int y0, int x1, int y1) {
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;
    int x = x0, y = y0;

    while (true) {
        if(is_blocked(x, y)) {
            return false;
        }
        if(x == x1 && y == y1) {
            break;
        }

        int e2 = 2 * err;

        if(e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if(e2 < dx) {
            err += dx;
            y += sy;
        }
    }
    return true;
}

/**
 * Very simple open list:
 * scan all nodes and pick the lowest f among opened and not closed
 * Slow, but easy to understand
 */
Node* get_best_open_node(void) {
    Node* best = NULL;
    for(int y = 0; y < MAX_ENVIRON_HEIGHT; y++) {
        for(int x = 0; x < MAX_ENVIRON_WIDTH; x++) {
            Node* n = &nodes[y][x];
            if(n->opened && !n->closed) {
                if(best == NULL || n->f < best->f) {
                    best = n;
                }
            }
        }
    }
    return best;
}

void init_nodes(int goal_x, int goal_y) {
    for(int y = 0; y < MAX_ENVIRON_HEIGHT; y++) {
        for(int x = 0; x < MAX_ENVIRON_WIDTH; x++) {
            nodes[y][x].x = x;
            nodes[y][x].y = y;
            nodes[y][x].g = INF;
            nodes[y][x].h = heuristic(x, y, goal_x, goal_y);
            nodes[y][x].f = INF;
            nodes[y][x].px = x;
            nodes[y][x].py = y;
            nodes[y][x].closed = false;
            nodes[y][x].opened = false;
        }
    }
}

void update_vertex(Node* current, Node* neighbor) {
    Node* parent = &nodes[current->py][current->px];

    if(line_of_sight(parent->x, parent->y, neighbor->x, neighbor->y)) {
        double new_g = parent->g + distance_between(parent->x, parent->y, neighbor->x, neighbor->y);
        if(new_g < neighbor->g) {
            neighbor->g = new_g;
            neighbor->px = parent->x;
            neighbor->py = parent->y;
            neighbor->f = neighbor->g + neighbor->h;
            neighbor->opened = true;
        }
    }
    else {
        double new_g = current->g + distance_between(current->x, current->y, neighbor->x, neighbor->y);
        if(new_g < neighbor->g) {
            neighbor->g = new_g;
            neighbor->px = current->x;
            neighbor->py = current->y;
            neighbor->f = neighbor->g + neighbor->h;
            neighbor->opened = true;
        }
    }
}

bool reconstruct_path(int start_x, int start_y, int goal_x, int goal_y, Point* path, int* path_length) {
    int max_path = MAX_ENVIRON_WIDTH * MAX_ENVIRON_HEIGHT;
    int count = 0;

    int cx = goal_x;
    int cy = goal_y;

    while(!(cx == start_x && cy == start_y)) {
        if(count >= max_path) {
            return false;
        }
        path[count++] = (Point){cx, cy};

        Node* n = &nodes[cy][cx];
        if(n->px == cx && n->py == cy) {
            return false;
        }

        int nx = n->px;
        int ny = n->py;
        cx = nx;
        cy = ny;
    }

    path[count++] = (Point){start_x, start_y};

    for(int i = 0; i < count / 2; i++) {
        Point temp = path[i];
        path[i] = path[count - 1 - i];
        path[count - 1 - i] = temp;
    }

    *path_length = count;
    return true;
}

bool theta_star(int start_x, int start_y, int goal_x, int goal_y, Point* path, int* path_length) {
    if(is_blocked(start_x, start_y) || is_blocked(goal_x, goal_y)) {
        return false;
    }

    init_nodes(goal_x, goal_y);

    Node* start = &nodes[start_y][start_x];
    start->g = 0.0;
    start->f = start->h;
    start->px = start_x;
    start->py = start_y;
    start->opened = true;

    while(true) {
        Node* current = get_best_open_node();
        if(current == NULL) {
            return false;
        }

        if(current->x == goal_x && current->y == goal_y) {
            return reconstruct_path(start_x, start_y, goal_x, goal_y, path, path_length);
        }

        current->closed = true;

        for(int i = 0; i < 8; i++) {
            int nx = current->x + dx[i];
            int ny = current->y + dy[i];

            if(!in_bounds(nx, ny) || is_blocked(nx, ny)) {
                continue;
            }

            Node* neighbor = &nodes[ny][nx];
            if(neighbor->closed) {
                continue;
            }

            if(!neighbor->opened) {
                neighbor->g = INF;
                neighbor->px = nx;
                neighbor->py = ny;
            }

            update_vertex(current, neighbor);
        }
    }
}

void print_grid_with_path(Point* path, int path_length, int start_x, int start_y, int goal_x, int goal_y) {
    char display[MAX_ENVIRON_HEIGHT][MAX_ENVIRON_WIDTH];

    for(int y = 0; y < MAX_ENVIRON_HEIGHT; y++) {
        for(int x = 0; x < MAX_ENVIRON_WIDTH; x++) {
            display[y][x] = (grid[y][x] == 0) ? '#' : '.';
        }
    }

    for(int i = 0; i < path_length; i++) {
        int x = path[i].x;
        int y = path[i].y;
        display[y][x] = '*';
    }

    display[start_y][start_x] = 'S';
    display[goal_y][goal_x] = 'G';

    for(int y = 0; y < MAX_ENVIRON_HEIGHT; y++) {
        for(int x = 0; x < MAX_ENVIRON_WIDTH; x++) {
            printf("%c ", display[y][x]);
        }
        printf("\n");
    }
}

int main(void) {
    int start_x = 0, start_y = 0;
    int goal_x = 9, goal_y = 9;

    Point path[MAX_ENVIRON_WIDTH * MAX_ENVIRON_HEIGHT];
    int path_length = 0;

    if(theta_star(start_x, start_y, goal_x, goal_y, path, &path_length)) {
        printf("Path found! Length: %d\n", path_length);
        print_grid_with_path(path, path_length, start_x, start_y, goal_x, goal_y);
    } else {
        printf("No path found.\n");
    }

    return 0;
}