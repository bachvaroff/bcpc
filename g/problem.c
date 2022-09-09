#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

enum state {
	INITIAL,
	RESTRICTION,
	POINT,
	FINAL
};

struct point {
	uint16_t x;
	uint16_t y;
};

struct Ys {
	uint16_t *Y;
	size_t N, size;
};

struct problem {
	struct point R;
	struct point alpha, omega;
	int concave;
	struct Ys X[USHRT_MAX];
	uint8_t *r;
	int N;
};

static void init_problem(struct problem *p) {
	(void)memset(p, 0, sizeof (*p));
	
	return;
}

static void cleanup(struct problem *p) {
	size_t j;
	
	for (j = 0u; j < USHRT_MAX; j++) {
		free(p->X[j].Y);
		p->X[j].Y = NULL;
		p->X[j].size = 0u;
		p->X[j].N = 0u;
	}
	
	free(p->r);
	
	return;
}

static void reset(struct problem *p) {
	size_t j;
	
	for (j = 0u; j < USHRT_MAX; j++)
		p->X[j].N = 0u;
	
	free(p->r);
	
	return;
}

static void parse_1(char *s, int *T) {
	char *p;
	
	p = strchr(s, (int)'\n');
	if (p) *p = 0;
	
	*T = atoi(s);
	
	return;
}

static void parse_2(char *s, uint16_t *x, uint16_t *y) {
	char *p;
	
	p = strchr(s, (int)'\n');
	if (p) *p = 0;
	
	assert((p = strchr(s, (int)' ')));
	
	*p++ = 0;
	
	*x = atoi(s);
	*y = atoi(p);
	
	return;
}

static void parse_3(char *s, uint16_t *x, uint16_t *y, int *N) {
	char *p, *t = s;
	
	if ((p = strchr(s, (int)'\n'))) *p = 0;
	
	assert((p = strchr(t, (int)' ')));
	*p++ = 0;
	*x = atoi(t);
	t = p;
	
	assert((p = strchr(t, (int)' ')));
	*p++ = 0;
	*y = atoi(t);
	t = p;
	
	assert((p = strchr(t, (int)' ')));
	*p++ = 0;
	*N = atoi(t);
	
	return;
}

static void do_realloc(struct problem *p, uint16_t x) {
	uint16_t *t;
	
	assert((t = realloc(p->X[x].Y, p->X[x].size + 1024lu * sizeof (uint16_t))));
	p->X[x].size += 1024lu * sizeof (uint16_t);
	p->X[x].Y = t;
	
	return;
}

static void update(struct problem *p, uint16_t x, uint16_t y) {
	if ((p->X[x].N + 1lu) >= (p->X[x].size / sizeof (uint16_t))) do_realloc(p, x);
	
	p->X[x].Y[p->X[x].N] = y;
	p->X[x].N++;
	
	return;
}

static int compare_u16(const void *_Y0, const void *_Y1) {
	const uint16_t *Y0 = _Y0, *Y1 = _Y1;
	
	if (*Y0 == *Y1) return 0;
	else if (*Y0 < *Y1) return -1;
	else return 1;
}

static void sort_Y(struct problem *p) {
	size_t x;
	
	for (x = 0u; x < USHRT_MAX; x++)
		if (p->X[x].Y) qsort(p->X[x].Y, p->X[x].N, sizeof (*p->X[x].Y), compare_u16);
	
	return;
}

static uint16_t *search_Y(const uint16_t *Y, const size_t N, const uint16_t _v) {
	uint16_t v = _v, *pv = &v;
	
	return bsearch(pv, Y, N, sizeof (*pv), compare_u16);
}

#define COORD(T, X, Y, XS) \
	((T)[(Y) * (XS) + (X)])

static uint8_t *init_restriction(struct problem *p) {
	uint16_t i, j;
	uint8_t *t;
	
	assert((t = (uint8_t *)malloc(p->R.x * p->R.y * sizeof (*t))));
	(void)memset(t, 0, p->R.x * p->R.y * sizeof (*t));
/*
	if (p->concave) {
		for (i = 1u; i < p->R.x; i++)
			for (j = 1u; j < p->R.y; j++)
				if (i < j) COORD(t, i, j, p->R.x) = 1;
	} else {
		for (i = 1u; i < p->R.x; i++)
			for (j = 1u; j < p->R.y; j++)
				if (j < i) COORD(t, i, j, p->R.x) = 1;
	}
*/
/*
	for (i = 1; i < p->R.x; i++)
		for (j = 1; j < p->R.y; j++)
			if ((i < j) || (i > j)) COORD(t, i, j, p->R.x) = 1;
*/
	for (i = 1u; i < p->R.x; i++)
		for (j = 1u; j < p->R.y; j++)
			COORD(t, i, j, p->R.x) = 1u;
	
	return t;
}

static int check(struct problem *p, struct point *three) {
	int32_t l, r;
	
	l = ((int32_t)three[2].y - (int32_t)three[1].y) * ((int32_t)three[1].x - (int32_t)three[0].x);
	r = ((int32_t)three[1].y - (int32_t)three[0].y) * ((int32_t)three[2].x - (int32_t)three[1].x);
		
	if (p->concave) {
		if (l < r) return 1;
	} else {
		if (l > r) return 1;
	}
	
	return 0;
}

static float dist(struct point *two) {
	float dx = (float)two[0].x - (float)two[1].x;
	float dy = (float)two[0].y - (float)two[1].y;
	
	return sqrtf(dx * dx + dy * dy);
}

static float max2(float max, float cd) {
	return (max > cd) ? max : cd;
}

static float recurse_points(struct problem *p, struct point *three, int valid, int depth) {
	struct point _three[3];
	uint32_t i, j;
	uint32_t nextx, nexty;
	float d, cd, maxd = -1.0f;
		
	if (!valid) {
		_three[2] = p->alpha;
		return recurse_points(p, _three, 1, depth + 1);
	} else {
		if ((three[2].x == p->omega.x) && (three[2].y == p->omega.y)) return 0.0f;
		for (i = 0u; i < p->R.y; i++) {
			nexty = (uint32_t)three[2].y + i;
			if (nexty > p->omega.y) continue;
			for (j = 0u; j < p->R.x; j++) {
				if (!COORD(p->r, j, i, p->R.x)) continue;
				nextx = (uint32_t)three[2].x + j;
				if (nextx > p->omega.x) continue;
				if (!p->X[nextx].Y) continue;
				if (search_Y(p->X[nextx].Y, p->X[nextx].N, (uint16_t)nexty)) {
					_three[0] = three[1];
					_three[1] = three[2];
					_three[2].x = nextx;
					_three[2].y = nexty;
					if (((valid >= 2) ? check(p, _three) : 1)) {
						cd = dist(_three + 1);
						d = recurse_points(p, _three, (valid >= 2) ? 3 : 2, depth + 1);
						if (d >= 0.0f) maxd = max2(maxd, cd + d);
					}
				}
			}
		}
	}
	
	return maxd;
}

static void solve(struct problem *p, FILE *output) {
	struct point three[3];
	float maxd;
		
	assert(p->X[p->alpha.x].Y);
	assert(p->X[p->omega.x].Y);
	
	sort_Y(p);
		
	assert((p->r = init_restriction(p)));
	
	maxd = recurse_points(p, three, 0, 0);
	
	fprintf(output, "%0.4f\n", maxd);
	
	return;
}

int main(void) {
	FILE *input, *output;
	char *line = NULL, *convex, *concave;
	size_t len;
	enum state s = INITIAL;
	struct problem *p;
	uint16_t x, y;
	int j;
	int T;
	
	assert((p = (struct problem *)malloc(sizeof (*p))));
	init_problem(p);
	
	input = stdin;
	output = stdout;
	
	while (getline(&line, &len, input) >= 0) {
		switch (s) {
		case INITIAL:
			parse_1(line, &T);
			s = RESTRICTION;
			break;
		case RESTRICTION:
			convex = strstr(line, "convex");
			concave = strstr(line, "concave");
			parse_3(line, &p->R.x, &p->R.y, &p->N);
			assert(p->R.x >= 3);
			assert(p->R.y >= 3);
			if (concave) p->concave = 1;
			else if (convex) p->concave = 0;
			else assert(0);
			s = POINT;
			j = 0;
			break;
		case POINT:
			parse_2(line, &x, &y);
			if (j == 0) {
				p->alpha.x = x;
				p->alpha.y = y;
			} else if (j == (p->N - 1)) {
				p->omega.x = x;
				p->omega.y = y;
				s = FINAL;
			}
			update(p, x, y);
			j++;
			break;
		default:
			break;
		}
		if (s == FINAL) {
			solve(p, output);
			if (--T) {
				reset(p);
				s = RESTRICTION;
			} else break;
		}
	}
	
	cleanup(p);
	free(p);
	free(line);
	
	return 0;
}
