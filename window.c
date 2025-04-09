/*!\file window.c
 * \brief GL4Dummies, exemple de physique
 * \author Farès BELHADJ, amsi@up8.edu
 * \date jannuary 20, 2025
 */

/* inclusion des entêtes de fonctions de création et de gestion de
 * fenêtres système ouvrant un contexte favorable à GL4dummies. Cette
 * partie est dépendante de la bibliothèque SDL2 */
#include <GL4D/gl4duw_SDL2.h>
#include <GL4D/gl4dm.h>
#include <GL4D/gl4dg.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NEIGHBOURS 64
#define HASH_SIZE 1024
#define CELL_SIZE 0.1f  //taille des cellules pour la grille spatiale

//pour la gravitée radiale
#define G_CONSTANT 6.67430e-3  // Constante gravitationnelle modifiée pour l'échelle de la simulation
#define MIN_DISTANCE 0.01f     // Distance minimale pour éviter les accélérations infinies
static int SPACE_MODE = 1;    // 0 = mode eau normal, 1 = mode espace gravitationnel
float vitesse = 1.0f;
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

//macro pour les expressions mathématiques complexes
#define POLY6 (315.0f / (64.0f * M_PI * powf(H, 9)))
#define SPIKY_GRAD (-45.0f / (M_PI * powf(H, 6)))
#define VISC_LAP (45.0f / (M_PI * powf(H, 6)))
#define SURFACE_TENSION 0.0728f

//pour les rectangles 3D
typedef struct {
    float x, y, z;
    float w, h, d;
    float angle_x; // Angle de rotation sur l'axe x
} rect3d_t;

//tableau dynamique global
static rect3d_t* _rects = NULL;  
static int _nb_rects   = 0;     
static int _max_rects  = 0;     


//struct pour la grille spatiale
typedef struct {
    int start_index;
    int count;
} spatial_cell_t;


//variables globales pour l'optimisation
static spatial_cell_t* _grid = NULL;
static int* _cell_indices = NULL;
static int* _particle_indices = NULL;


typedef struct vec3d_t vec3d_t;
typedef struct mobile_t mobile_t;

struct vec3d_t
{
	GLfloat x, y, z;
};

//struct mobile_t
//{
//	vec3d_t p, v;
//	GLfloat r;
//	GLfloat color[4];
//};

//temps de simulation
static float TIME_SCALE = 10.0f;  //vitesse de simulation


static void init(void);
static void draw(void);
static void quit(void);

static void mobile_init(int n);
static void mobile_simu(void);
static void mobile_draw(void);
static void mobile_quit(void);

static void mobile_simu_with_mode_selection(void);

/* on créé une variable pour stocker l'identifiant du programme GPU */
GLuint _pId = 0;

GLuint _quad = 0;

/* gravité */
//static GLfloat _ig = 9.81f / 2.0f;
static vec3d_t _g = {0.0f, 0.0f, 0.0f}; // Modification ici: définir la gravité vers le bas à -9.81f
static const GLfloat e = 0.5f; //8.0f / 9.0f;

/* simulation d'eau de jsp qui */
// Ajouter ces paramètres SPH
static const float REST_DENSITY = 300.0f;  // Densité au repos du fluide
static const float GAS_CONSTANT = 2000.0f;  // Constante des gaz parfaits
static const float VISCOSITY = 10.0f;      // Viscosité du fluide
static const float MASS = 1.0f;             // Masse d'une particule
static const float H = 0.11f;                // Rayon de lissage (smoothing radius)
static const float H2 = 0.0075f;              // H²
//static const float POLY6 = 315.0f / (64.0f * M_PI * powf(H, 9));
//static const float SPIKY_GRAD = -45.0f / (M_PI * powf(H, 6));
//static const float VISC_LAP = 45.0f / (M_PI * powf(H, 6));
//static const float SURFACE_TENSION = 0.0728f;

// Variables fluides à ajouter à la structure mobile_t
struct mobile_t {
    vec3d_t p, v;
    GLfloat r;
    GLfloat color[4];
    
    // Variables SPH
    float density;
    float pressure;
    vec3d_t force;
    int cell_id;
};


/* tous les mobiles de ma scène */
static mobile_t *_mobiles = NULL;
static int _nb_mobiles = 0;

static int _ww = 600, _wh = 600;

/*!\brief créé la fenêtre, un screen 2D effacé en noir et lance une
 *  boucle infinie.*/
int main(int argc, char **argv){
	/* tentative de création d'une fenêtre pour GL4Dummies */
	if (!gl4duwCreateWindow(argc, argv,				   /* args du programme */
							"GL4Dummies' Hello World", /* titre */
							10, 10, _ww, _wh,		   /* x,y, largeur, heuteur */
							GL4DW_SHOWN) /* état visible */){
		/* ici si échec de la création souvent lié à un problème d'absence
		 * de contexte graphique ou d'impossibilité d'ouverture d'un
		 * contexte OpenGL (au moins 3.2) */
		return 1;
	}
	/* appeler init pour initialiser des paramètres GL et GL4D */
	init();
	/* placer quit comme fonction à appeler au moment du exit */
	atexit(quit);
	/* placer mobile_simu comme fonction à appeler à idle (simulation) */
	gl4duwIdleFunc(mobile_simu);
    //gl4duwIdleFunc(mobile_simu_with_mode_selection);
	/* placer draw comme fonction à appeler pour dessiner chaque frame */
	gl4duwDisplayFunc(draw);
	/* boucle infinie pour éviter que le programme ne s'arrête et ferme
	 * la fenêtre immédiatement */
	gl4duwMainLoop();
	return 0;
}

// Fonctions auxiliaires pour SPH
float kernel_poly6(float r2) {
    if (r2 > H2) return 0.0f;
    float temp = H2 - r2;
    return POLY6 * temp * temp * temp;
}

float kernel_spiky_gradient(float r) {
    if (r > H) return 0.0f;
    float temp = H - r;
    return SPIKY_GRAD * temp * temp;
}

float kernel_viscosity_laplacian(float r) {
    if (r > H) return 0.0f;
    return VISC_LAP * (H - r);
}
float kernel_viscosity_improved(float r, float h) {
    if (r >= h) return 0.0f;
    float q = r / h;
    return (h / r) * (1.0f - q);
}

//les rectangles
void rect_init_list(int capacity) {
    if (_rects) free(_rects);
    _max_rects = capacity;
    _rects = (rect3d_t*)malloc(_max_rects * sizeof(rect3d_t));
    _nb_rects = 0;
}
// Ajoute un rectangle à la liste
void rect_add(float x, float y, float z, float w, float h, float d, float angle_x) {
    if (_nb_rects >= _max_rects) return; // ou agrandir dynamiquement
    _rects[_nb_rects].x = x;
    _rects[_nb_rects].y = y;
    _rects[_nb_rects].z = z;
    _rects[_nb_rects].w = w;
    _rects[_nb_rects].h = h;
    _rects[_nb_rects].d = d;
    _rects[_nb_rects].angle_x = angle_x; // Ajout de l'angle de rotation sur l'axe x
    _nb_rects++;
}

// Test collision (2D) d'une particule (mobile_t) avec un rectangle avec rotation
static void collide_with_rotated_rect(mobile_t* m, const rect3d_t* r, float e) {
    // Translate the particle's position to the rectangle's local space
    float localX = m->p.x - r->x;
    float localY = m->p.y - r->y;

    // Apply rotation
    float cosAngle = cosf(r->angle_x);
    float sinAngle = sinf(r->angle_x);
    float rotatedX = cosAngle * localX + sinAngle * localY;
    float rotatedY = -sinAngle * localX + cosAngle * localY;

    // Check if the rotated coordinate is within the rectangle's bounds
    if (rotatedX >= 0.0f && rotatedX <= r->w &&
        rotatedY >= 0.0f && rotatedY <= r->h) {
        
        // Correct the particle's position to push it outside the rectangle
        if (rotatedY < r->h / 2.0f) {
            rotatedY = -m->r; // Push below
        } else {
            rotatedY = r->h + m->r; // Push above
        }

        // Apply inverse rotation to return to global space
        float correctedX = cosAngle * rotatedX - sinAngle * rotatedY;
        float correctedY = sinAngle * rotatedX + cosAngle * rotatedY;

        // Update the particle's position and velocity
        m->p.x = correctedX + r->x;
        m->p.y = correctedY + r->y;
        m->v.y = -m->v.y * e;
    }
}

// Dessine tous les rectangles en utilisant le shader program
void rect_draw_all(void) {
    GLfloat *rect_data = malloc(4 * _nb_rects * sizeof *rect_data); //4 pour x,y,w,h
    assert(rect_data);
    
    for (int i = 0; i < _nb_rects; ++i) {
        rect_data[4 * i + 0] = _rects[i].x;
        rect_data[4 * i + 1] = _rects[i].y;
        rect_data[4 * i + 2] = _rects[i].w;
        rect_data[4 * i + 3] = _rects[i].h;
    }
    // Set up an attribute for angles if needed
    glUseProgram(_pId);
    //allocation séparée pour les angles
    float *angles = malloc(_nb_rects * sizeof(float));
    assert(angles);
    for (int i = 0; i < _nb_rects; ++i) {
        angles[i] = _rects[i].angle_x;
    }
    glUniform1fv(glGetUniformLocation(_pId, "rect_angles"), _nb_rects, angles);
    glUseProgram(_pId);
    glUniform4fv(glGetUniformLocation(_pId, "rectangles"), _nb_rects, rect_data);
    glUniform1i(glGetUniformLocation(_pId, "nb_rects"), _nb_rects);

    //set la couleur en blanc
    glUniform4f(glGetUniformLocation(_pId, "rect_color"), 1.0f, 1.0f, 1.0f, 1.0f);
    for (int i = 0; i < _nb_mobiles; ++i) {
        for (int j = 0; j < _nb_rects; ++j) {
            collide_with_rotated_rect(&_mobiles[i], &_rects[j], e);
        }
    }
    gl4dgDraw(_quad);
    glUseProgram(0);

    free(rect_data);
    free(angles);
}


// Appeler cette fonction dans mobile_simu après mise à jour des particules
void rect_collide_all(mobile_t* mobiles, int nb, float e) {
    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < _nb_rects; j++) {
            collide_with_rotated_rect(&mobiles[i], &_rects[j], e);
        }
    }
}

// Libère la liste
void rect_cleanup(void) {
    if (_rects) {
        free(_rects);
        _rects = NULL;
    }
    _nb_rects = 0;
    _max_rects = 0;
}


// Fonction pour construire la grille spatiale
void build_spatial_grid() {
    // Initialiser la grille
    if (!_grid) {
        _grid = (spatial_cell_t*)malloc(HASH_SIZE * sizeof(spatial_cell_t));
        _cell_indices = (int*)malloc(_nb_mobiles * sizeof(int));
        _particle_indices = (int*)malloc(_nb_mobiles * sizeof(int));
    }
    
    // Réinitialiser les cellules
    for (int i = 0; i < HASH_SIZE; i++) {
        _grid[i].start_index = -1;
        _grid[i].count = 0;
    }
    
    // Attribuer les particules aux cellules
    for (int i = 0; i < _nb_mobiles; i++) {
        int cellX = (int)(((_mobiles[i].p.x + 1.0f) / 2.0f) / CELL_SIZE);
        int cellY = (int)(((_mobiles[i].p.y + 1.0f) / 2.0f) / CELL_SIZE);
        
        int cell_id = (cellY * 32 + cellX) % HASH_SIZE; // Hachage simple
        _mobiles[i].cell_id = cell_id;
        
        // Compter combien de particules par cellule
        _grid[cell_id].count++;
    }
    
    // Calculer les indices de départ
    int start = 0;
    for (int i = 0; i < HASH_SIZE; i++) {
        _grid[i].start_index = start;
        start += _grid[i].count;
        _grid[i].count = 0; // Réinitialiser pour la prochaine étape
    }
    
    // Remplir les indices
    for (int i = 0; i < _nb_mobiles; i++) {
        int cell_id = _mobiles[i].cell_id;
        int index = _grid[cell_id].start_index + _grid[cell_id].count;
        _cell_indices[index] = i;
        _grid[cell_id].count++;
    }
}

// Fonction pour calculer les forces gravitationnelles entre particules
void compute_gravity_forces() {
    //reset les forces
    for (int i = 0; i < _nb_mobiles; i++) {
        _mobiles[i].force.x = 0.0f;
        _mobiles[i].force.y = 0.0f;
        _mobiles[i].force.z = 0.0f;
    }
    
    //calcule des forces gravitationnelles entre chaque paire de particules
    for (int i = 0; i < _nb_mobiles; i++) {
        for (int j = i + 1; j < _nb_mobiles; j++) {
            float dx = _mobiles[j].p.x - _mobiles[i].p.x;
            float dy = _mobiles[j].p.y - _mobiles[i].p.y;
            //float dz = _mobiles[j].p.z - _mobiles[i].p.z;
            
            float r2 = dx*dx + dy*dy;// + dz*dz;
            float r = sqrtf(r2);
            
            // Éviter division par zéro et forces trop grandes entre particules proches
            if (r < MIN_DISTANCE) r = MIN_DISTANCE;
            
            //force gravitationnelle: F = G * m1 * m2 / r^2
            //toutes les particules ont la même masse (MASS)
            float force = G_CONSTANT * MASS * MASS / r2;
            
            //direction de la force (vecteur unitaire)
            float fx = force * dx / r;
            float fy = force * dy / r;
            //float fz = force * dz / r;
            
            //la force aux deux particules (action-réaction)
            _mobiles[i].force.x += fx;
            _mobiles[i].force.y += fy;
            //_mobiles[i].force.z += fz;
            
            _mobiles[j].force.x -= fx;
            _mobiles[j].force.y -= fy;
            //_mobiles[j].force.z -= fz;

            // Réduire la vitesse des particules pour simuler un amortissement
            _mobiles[i].v.x *= 0.9f; // Réduction de 80% par itération
            _mobiles[i].v.y *= 0.9f;
            //_mobiles[i].v.z *= 0.9f;
            _mobiles[j].v.x *= 0.9f;
            _mobiles[j].v.y *= 0.9f;
            //_mobiles[j].v.z *= 0.9f;
            _mobiles[i].force.x *= 0.1f;
            _mobiles[i].force.y *= 0.1f;
        }
    }
}


// Gérer les collisions physiques entre particules
void handle_particle_collisions(float dt) {
    for (int i = 0; i < _nb_mobiles; ++i) {
        for (int j = i + 1; j < _nb_mobiles; ++j) {
            float dx = _mobiles[j].p.x - _mobiles[i].p.x;
            float dy = _mobiles[j].p.y - _mobiles[i].p.y;
            float dz = _mobiles[j].p.z - _mobiles[i].p.z;
            
            float dist2 = dx*dx + dy*dy + dz*dz;
            float min_dist = _mobiles[i].r + _mobiles[j].r;
            
            if (dist2 < min_dist * min_dist) {
                float dist = sqrtf(dist2);
                if (dist < 0.0001f) dist = 0.0001f; // Éviter division par zéro
                
                // Vecteurs unitaires
                float nx = dx / dist;
                float ny = dy / dist;
                float nz = dz / dist;
                
                // Correction de position pour éviter le chevauchement
                float overlap = 0.5f * (min_dist - dist);
                _mobiles[i].p.x -= nx * overlap;
                _mobiles[i].p.y -= ny * overlap;
                _mobiles[i].p.z -= nz * overlap;
                _mobiles[j].p.x += nx * overlap;
                _mobiles[j].p.y += ny * overlap;
                _mobiles[j].p.z += nz * overlap;
                
                // Calcul de l'impulsion pour le rebond
                float vx = _mobiles[j].v.x - _mobiles[i].v.x;
                float vy = _mobiles[j].v.y - _mobiles[i].v.y;
                float vz = _mobiles[j].v.z - _mobiles[i].v.z;
                
                float dot = vx*nx + vy*ny + vz*nz;
                if (dot > 0.0f) continue; // Les particules s'éloignent déjà
                
                float impulse = -(1.0f + e) * dot;
                impulse /= 2.0f; // Masses égales
                
                _mobiles[i].v.x -= impulse * nx;
                _mobiles[i].v.y -= impulse * ny;
                _mobiles[i].v.z -= impulse * nz;
                
                _mobiles[j].v.x += impulse * nx;
                _mobiles[j].v.y += impulse * ny;
                _mobiles[j].v.z += impulse * nz;
            }
        }
    }
}


// Ajouter une nouvelle fonction de simulation pour le mode spatial
void space_simulation() {
    static double t0 = 0;
    double t = gl4dGetElapsedTime() / 1000.0, dt = (t - t0) * 30.0;
    t0 = t;
    
    if (dt > 0.03f) dt = 0.03f; // Limiter le pas de temps
    
    // Calculer les forces gravitationnelles entre particules
    compute_gravity_forces();
    
    // Mise à jour des positions et vitesses
    for (int i = 0; i < _nb_mobiles; ++i) {
        // Intégration d'Euler
        _mobiles[i].v.x += _mobiles[i].force.x * dt / MASS;
        _mobiles[i].v.y += _mobiles[i].force.y * dt / MASS;
        _mobiles[i].v.z += _mobiles[i].force.z * dt / MASS;
        
        _mobiles[i].p.x += _mobiles[i].v.x * dt;
        _mobiles[i].p.y += _mobiles[i].v.y * dt;
        _mobiles[i].p.z += _mobiles[i].v.z * dt;
        
        // Limiter la vitesse maximale pour la stabilité
        float speed = sqrtf(_mobiles[i].v.x * _mobiles[i].v.x + 
                            _mobiles[i].v.y * _mobiles[i].v.y +
                            _mobiles[i].v.z * _mobiles[i].v.z);
        const float max_speed = 1.0f;
        if (speed > max_speed) {
            float ratio = max_speed / speed;
            _mobiles[i].v.x *= ratio;
            _mobiles[i].v.y *= ratio;
            _mobiles[i].v.z *= ratio;
        }
        
        // Rebond sur les bords (ou espace toroïdal)
        if (_mobiles[i].p.x < -1.0f) {
            _mobiles[i].p.x = -1.0f;
            _mobiles[i].v.x = -_mobiles[i].v.x * e;
        }
        if (_mobiles[i].p.x > 1.0f) {
            _mobiles[i].p.x = 1.0f;
            _mobiles[i].v.x = -_mobiles[i].v.x * e;
        }
        if (_mobiles[i].p.y < -1.0f) {
            _mobiles[i].p.y = -1.0f;
            _mobiles[i].v.y = -_mobiles[i].v.y * e;
        }
        if (_mobiles[i].p.y > 1.0f) {
            _mobiles[i].p.y = 1.0f;
            _mobiles[i].v.y = -_mobiles[i].v.y * e;
        }
        
        // Gestion du Z pour une éventuelle visualisation 3D future
        if (_mobiles[i].p.z < -1.0f) {
            _mobiles[i].p.z = -1.0f;
            _mobiles[i].v.z = -_mobiles[i].v.z * e;
        }
        if (_mobiles[i].p.z > 1.0f) {
            _mobiles[i].p.z = 1.0f;
            _mobiles[i].v.z = -_mobiles[i].v.z * e;
        }
        
        // Coloration basée sur la vitesse pour visualiser l'énergie
        float energy = speed / max_speed;
        _mobiles[i].color[0] = energy;
        _mobiles[i].color[1] = 0.4f;
        _mobiles[i].color[2] = 1.0f - energy * 0.5f;
    }
    
    // Gérer les collisions entre particules
    handle_particle_collisions(dt);
    
    // Gérer les collisions avec les obstacles
    rect_collide_all(_mobiles, _nb_mobiles, e);
}

// Initialisation pour mode spatial
void space_init(int n) {
    assert(_mobiles == NULL);
    _nb_mobiles = n;
    _mobiles = malloc(_nb_mobiles * sizeof *_mobiles);
    assert(_mobiles);
    
    //configuration initiale en forme de disque/anneau
    float center_x = 0.0f;
    float center_y = 0.0f;
    float inner_radius = 0.1f;
    float outer_radius = 0.7f;
    float orbital_velocity = 0.4f;  // Vitesse orbitale de base

    for (int i = 0; i < n; i++) {
        // Position aléatoire dans un anneau
        float angle = gl4dmSURand() * 2.0f * M_PI;
        float radius = inner_radius + gl4dmSURand() * (outer_radius - inner_radius);
        
        _mobiles[i].p.x = center_x + radius * cosf(angle);
        _mobiles[i].p.y = center_y + radius * sinf(angle);
        _mobiles[i].p.z = (gl4dmSURand() - 0.5f) * 0.2f;  // Petite variation en Z
        
        // Vitesse orbitale (perpendiculaire au rayon)
        float orbital_speed = orbital_velocity / sqrtf(radius);  // Loi de Kepler simplifiée
        _mobiles[i].v.x = -sinf(angle) * orbital_speed * (0.8f + 0.4f * gl4dmSURand());
        _mobiles[i].v.y = cosf(angle) * orbital_speed * (0.8f + 0.4f * gl4dmSURand());
        _mobiles[i].v.z = 0.0f;
        
        // Taille et autres propriétés
        _mobiles[i].r = 0.01f + 0.01f * gl4dmSURand();  // Tailles variables
        _mobiles[i].density = REST_DENSITY;
        _mobiles[i].pressure = 0.0f;
        _mobiles[i].force.x = _mobiles[i].force.y = _mobiles[i].force.z = 0.0f;
        
        // Couleurs basées sur la distance au centre (comme les anneaux planétaires)
        float color_factor = (radius - inner_radius) / (outer_radius - inner_radius);
        _mobiles[i].color[0] = 0.8f - 0.4f * color_factor;
        _mobiles[i].color[1] = 0.4f + 0.4f * color_factor;
        _mobiles[i].color[2] = 0.6f + 0.4f * color_factor;
        _mobiles[i].color[3] = 1.0f;
    }
    
    // Ajouter un corps central plus massif (ajustez l'indice si nécessaire)
    if (n > 0) {
        _mobiles[0].p.x = center_x;
        _mobiles[0].p.y = center_y;
        _mobiles[0].p.z = 0.0f;
        _mobiles[0].v.x = _mobiles[0].v.y = _mobiles[0].v.z = 0.0f;
        _mobiles[0].r = 0.05f;  // Corps central plus grand
        _mobiles[0].color[0] = 1.0f;
        _mobiles[0].color[1] = 0.9f;
        _mobiles[0].color[2] = 0.4f;
    }
}

// Fonction principale SPH
void compute_sph_forces() {
    // Construire la grille spatiale
    build_spatial_grid();
    // Calculer les densités et pressions
    for (int i = 0; i < _nb_mobiles; i++) {
        _mobiles[i].density = 0.0f;
        
        // Auto-contribution à la densité
        _mobiles[i].density += MASS * kernel_poly6(0.0f);
        
        // Contribution des voisins
        //int cell_id = _mobiles[i].cell_id;
        
        // Parcourir les 9 cellules voisines (en 2D)
        for (int offsetY = -1; offsetY <= 1; offsetY++) {
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                int cellX = (int)(((_mobiles[i].p.x + 1.0f) / 2.0f) / CELL_SIZE) + offsetX;
                int cellY = (int)(((_mobiles[i].p.y + 1.0f) / 2.0f) / CELL_SIZE) + offsetY;
                
                if (cellX < 0 || cellY < 0 || cellX >= 32 || cellY >= 32) continue;
                
                int neighbor_cell_id = (cellY * 32 + cellX) % HASH_SIZE;
                
                if (_grid[neighbor_cell_id].start_index == -1) continue;
                
                // Parcourir les particules de cette cellule
                for (int k = 0; k < _grid[neighbor_cell_id].count; k++) {
                    int j = _cell_indices[_grid[neighbor_cell_id].start_index + k];
                    if (i == j) continue;
                    
                    float dx = _mobiles[j].p.x - _mobiles[i].p.x;
                    float dy = _mobiles[j].p.y - _mobiles[i].p.y;
                    float r2 = dx*dx + dy*dy;
                    
                    if (r2 < H2) {
                        _mobiles[i].density += MASS * kernel_poly6(r2);
                    }
                }
            }
        }
        
        // Calculer la pression avec l'équation d'état
        //_mobiles[i].pressure = GAS_CONSTANT * (_mobiles[i].density - REST_DENSITY);
		float density_ratio = _mobiles[i].density / REST_DENSITY;
		if (density_ratio > 1.0f) {
		    // Pression positive (répulsive) si densité > REST_DENSITY
		    _mobiles[i].pressure = GAS_CONSTANT * (powf(density_ratio, 4) - 1.0f);
		} else {
		    // Pression négative (attractive) si densité < REST_DENSITY, mais plus faible
		    //ici on gere le mode espace
            if(SPACE_MODE){
                _mobiles[i].pressure = GAS_CONSTANT * 200.0f * (density_ratio - 1.0f);    
                if (density_ratio < 1.0f) {
                    // Attraction même à faible densité
                    _mobiles[i].pressure += GAS_CONSTANT * 10.0f * (1.0f - density_ratio);
                }
            }else{
                _mobiles[i].pressure = GAS_CONSTANT * 2.0f * (density_ratio - 1.0f);
            }
		}
		//float density_ratio = _mobiles[i].density / REST_DENSITY;
        //_mobiles[i].pressure = GAS_CONSTANT * (powf(density_ratio, 7) - 1.0f);
        
        // Limiter les pressions extrêmes pour éviter l'instabilité
        //if (_mobiles[i].pressure > GAS_CONSTANT * 10.0f)
        //    _mobiles[i].pressure = GAS_CONSTANT * 10.0f;
        //if (_mobiles[i].pressure < -GAS_CONSTANT)
        //    _mobiles[i].pressure = -GAS_CONSTANT;
        
    }
    
    // Calculer les forces
    for (int i = 0; i < _nb_mobiles; i++) {
        
        //int cell_id = _mobiles[i].cell_id;
        
        // Parcourir les 9 cellules voisines
        for (int offsetY = -1; offsetY <= 1; offsetY++) {
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                int cellX = (int)(((_mobiles[i].p.x + 1.0f) / 2.0f) / CELL_SIZE) + offsetX;
                int cellY = (int)(((_mobiles[i].p.y + 1.0f) / 2.0f) / CELL_SIZE) + offsetY;
                
                if (cellX < 0 || cellY < 0 || cellX >= 32 || cellY >= 32) continue;
                
                int neighbor_cell_id = (cellY * 32 + cellX) % HASH_SIZE;
                
                if (_grid[neighbor_cell_id].start_index == -1) continue;
                
                // Parcourir les particules de cette cellule
                for (int k = 0; k < _grid[neighbor_cell_id].count; k++) {
                    int j = _cell_indices[_grid[neighbor_cell_id].start_index + k];
                    if (i == j) continue;
                    
                    float dx = _mobiles[j].p.x - _mobiles[i].p.x;
                    float dy = _mobiles[j].p.y - _mobiles[i].p.y;
                    float r2 = dx*dx + dy*dy;
                    float r = sqrtf(r2);
                    
                    if (r < H && r > 0.0001f) {
                        // Force de pression
                        float pressure_factor = -MASS * (_mobiles[i].pressure + _mobiles[j].pressure) / 
                                               (2.0f * _mobiles[j].density) * kernel_spiky_gradient(r);
                        //printf("Pressure factor: %f\n", pressure_factor);
                        _mobiles[i].force.x += pressure_factor * dx / r;
                        _mobiles[i].force.y += pressure_factor * dy / r;
                        
                        // Force de viscosité
                        float visc_factor = VISCOSITY * MASS * 
                                         (_mobiles[j].v.x - _mobiles[i].v.x) * 
                                         //kernel_viscosity_laplacian(r) / _mobiles[j].density;
										 kernel_viscosity_improved(r, H);
                        
                        _mobiles[i].force.x += visc_factor*0.1f;
                        
                        visc_factor = VISCOSITY * MASS * 
                                   (_mobiles[j].v.y - _mobiles[i].v.y) * 
                                   //kernel_viscosity_laplacian(r) / _mobiles[j].density;
								   kernel_viscosity_improved(r, H);
                        
                        _mobiles[i].force.y += visc_factor;
						float min_distance = _mobiles[i].r + _mobiles[j].r;
						// Force de répulsion à courte distance (deux-couches)
						if (r < min_distance * 1.5f) {
						    // Première couche - très forte répulsion si presque en contact
						    if (r < min_distance * 1.1f) {
						        float repulsion_strength = 50.0f * (min_distance * 1.1f - r) / min_distance;
						        _mobiles[i].force.x += repulsion_strength * dx / r;
						        _mobiles[i].force.y += repulsion_strength * dy / r;
						    } 
						    // Deuxième couche - répulsion plus douce
						    else {
						        float repulsion_strength = 5.0f * (min_distance * 1.5f - r) / min_distance;
						        _mobiles[i].force.x += repulsion_strength * dx / r;
						        _mobiles[i].force.y += repulsion_strength * dy / r;
						    }
						}
						min_distance = _mobiles[i].r + _mobiles[j].r;
						if (r < min_distance) {
							float repulsion_strength = 10.0f * (min_distance - r) / min_distance;
							_mobiles[i].force.x += repulsion_strength * dx / r;
							_mobiles[i].force.y += repulsion_strength * dy / r;
						}
						
						// Limiter la densité
						if (_mobiles[i].density > REST_DENSITY * 2.0f) {
							_mobiles[i].v.x *= 0.5f;
							_mobiles[i].v.y *= 0.5f;
						}
						if (r < min_distance * 1.5f) {
						    float repulsion_strength = 2.0f * (min_distance * 1.5f - r) / min_distance;
						    _mobiles[i].force.x -= repulsion_strength * dx / r;
						    _mobiles[i].force.y -= repulsion_strength * dy / r;
						}
                    }
                }
            }
        }
		//limiter leur magnitude
        float force_magnitude = sqrtf(_mobiles[i].force.x * _mobiles[i].force.x + 
                                     _mobiles[i].force.y * _mobiles[i].force.y);
        const float max_force = 1000.0f;
        
        if (force_magnitude > max_force) {
            float scale = max_force / force_magnitude;
            _mobiles[i].force.x *= scale;
            _mobiles[i].force.y *= scale;
        }
        //_mobiles[i].force.x *= 0.1f;
        //_mobiles[i].force.y *= 0.1f;
		_mobiles[i].force.x += _g.x * _mobiles[i].density * 0.70f;
		_mobiles[i].force.y += _g.y * _mobiles[i].density * 0.70f;
		//_mobiles[i].force.x += _g.x;  // Indépendant de la densité
		//_mobiles[i].force.y += _g.y;  // Indépendant de la densité
        //printf("Force: (%f, %f)\n", _mobiles[i].force.x, _mobiles[i].force.y);
    }
}
/*fin fonction de ouf*/

/* initialise des paramètres GL et GL4D */
void init(void){
	_quad = gl4dgGenQuadf();
	/* activer la synchronisation verticale */
	SDL_GL_SetSwapInterval(1);
	/* set la couleur d'effacement OpenGL */
    if(!SPACE_MODE){
        TIME_SCALE = 10.0f;
        vitesse = 10.0f;
        _g.y = -9.81f;
    }else{
        TIME_SCALE = 10.0f;
        vitesse = 100.0f;
        _g.y = 0.0f;
        _g.x = 0.01f;
    }
    glClearColor(0.0f, 0.0f, 0.5f, 1.0f);
    /* créer un programme GPU pour OpenGL (en GL4D) */
    _pId = gl4duCreateProgram("<vs>shaders/identity.vs", "<fs>shaders/calculs.fs", NULL);
    mobile_init(1023);
    //les rectangles
    rect_init_list(3); //la liste de rectangles
    rect_add(0.50f, -0.0f, 0.0f, 1.10f, 0.10f, 0.10f, -3.01f); // Rectangle 1
    rect_add(0.4f, -0.1f, 0.0f, 0.10f, 1.10f, 0.10f, 0.0f); // Rectangle 2
    rect_add(-0.1f, 0.5f, 0.0f, 1.0f, 0.10f, 0.10f, 3.0f); // Rectangle 3
    //}else{
    //    glClearColor(0.5f, 0.0f, 0.0f, 1.0f);
    //    /* créer un programme GPU pour OpenGL (en GL4D) */
    //    _pId = gl4duCreateProgram("<vs>shaders/identity.vs", "<fs>shaders/calculs.fs", NULL);
    //    space_init(1023);
    //}
}


void draw(void){
	/* effacer le buffer de couleur (image) et le buffer de profondeur d'OpenGL */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* utiliser le programme GPU "_pId" */
	glUseProgram(_pId);
	/* binder (mettre au premier plan, "en courante" ou "en active") la
	   matrice view */
	mobile_draw();
    rect_draw_all();
    /* n'utiliser aucun programme GPU (pas nécessaire) */
    glUseProgram(0);
}

/* appelée lors du exit */
void quit(void){
	mobile_quit();
	/* nettoyer (libérer) tout objet créé avec GL4D */
	gl4duClean(GL4DU_ALL);
}

void mobile_init(int n){
	assert(_mobiles == NULL);
    _nb_mobiles = n;
    _mobiles = malloc(_nb_mobiles * sizeof *_mobiles);
    assert(_mobiles);
    // Créer une forme de carré pour l'initialisation
    float start_x = -0.5f;
    float start_y = 0.5f;
    float width = 0.6f;
    float height = 0.6f;
    float spacing = 0.05f;

    int index = 0;
    for (float y = start_y; y <= start_y + height && index < n; y += spacing) {
        for (float x = start_x; x <= start_x + width && index < n; x += spacing) {
            _mobiles[index].p.x = x;
            _mobiles[index].p.y = y;
            _mobiles[index].p.z = 0.0f;
            _mobiles[index].v.x = 0.0f;
            _mobiles[index].v.y = 0.0f;
            _mobiles[index].v.z = 0.0f;
            _mobiles[index].r = 0.01f;
            _mobiles[index].density = REST_DENSITY;
            _mobiles[index].pressure = 0.0f;
            _mobiles[index].force.x = 0.0f;
            _mobiles[index].force.y = 0.0f;
            _mobiles[index].force.z = 0.0f;
            _mobiles[index].color[0] = 0.0f;
            _mobiles[index].color[1] = 0.4f;
            _mobiles[index].color[2] = 0.8f;
            _mobiles[index].color[3] = 1.0f;
            index++;
        }
    }

    /*
    // Créer une forme de "goutte d'eau" pour l'initialisation
    float center_x = -0.5f;
    float center_y = 0.7f;
    float radius = 0.3f;
    float spacing = 0.05f;
    
    int index = 0;
    for (float y = center_y - radius; y <= center_y + radius && index < n; y += spacing) {
        for (float x = center_x - radius; x <= center_x + radius && index < n; x += spacing) {
            float dx = x - center_x;
            float dy = y - center_y;
            float dist = sqrtf(dx*dx + dy*dy);
            
            if (dist <= radius) {
                _mobiles[index].p.x = x;
                _mobiles[index].p.y = y;
                _mobiles[index].p.z = 0.0f;
                _mobiles[index].v.x = 0.0f;
                _mobiles[index].v.y = 0.0f;
                _mobiles[index].v.z = 0.0f;
                _mobiles[index].r = 0.02f;
                _mobiles[index].density = REST_DENSITY;
                _mobiles[index].pressure = 0.0f;
                _mobiles[index].force.x = 0.0f;
                _mobiles[index].force.y = 0.0f;
                _mobiles[index].force.z = 0.0f;
                _mobiles[index].color[0] = 0.0f;
                _mobiles[index].color[1] = 0.4f;
                _mobiles[index].color[2] = 0.8f;
                _mobiles[index].color[3] = 1.0f;
                index++;
            }
        }
    }
    */
    
    // Si on n'a pas assez de particules, remplir le reste
    for (int i = index; i < n; i++) {
        _mobiles[i].p.x = gl4dmSURand() * 0.5f - 0.5f;
        _mobiles[i].p.y = gl4dmSURand() * 0.5f + 0.5f;
        _mobiles[i].p.z = 0.0f;
        _mobiles[i].v.x = 0.0f;
        _mobiles[i].v.y = 0.0f;
        _mobiles[i].v.z = 0.0f;
        _mobiles[i].r = 0.01f;
        _mobiles[i].density = REST_DENSITY;
        _mobiles[i].pressure = 0.0f;
        _mobiles[i].force.x = 0.0f;
        _mobiles[i].force.y = 0.0f;
        _mobiles[i].force.z = 0.0f;
        _mobiles[i].color[0] = 0.0f;
        _mobiles[i].color[1] = 0.4f;
        _mobiles[i].color[2] = 0.8f;
        _mobiles[i].color[3] = 1.0f;
    }
}



void mobile_simu(void){
	static double t0 = 0;
	double t = gl4dGetElapsedTime() / 1000.0, dt = (t - t0) * 30.0;
	t0 = t;

	if (dt > 0.03f) dt = 0.03f; // Limiter le pas de temps à 30 ms
    // Calculer les forces SPH
    compute_sph_forces();
    for (int i = 0; i < _nb_mobiles; ++i) {
        _mobiles[i].force.x *= dt*TIME_SCALE;
        _mobiles[i].force.y *= dt*TIME_SCALE;
        _mobiles[i].force.z *= dt*TIME_SCALE;
    }
	for (int i = 0; i < _nb_mobiles; ++i) {
        // Intégration explicite d'Euler
        float accel_x = _mobiles[i].force.x / _mobiles[i].density;
        float accel_y = _mobiles[i].force.y / _mobiles[i].density;
        
        _mobiles[i].v.x += accel_x * dt * vitesse;
        _mobiles[i].v.y += accel_y * dt * vitesse;//TODO aller plus vite dans l'espace
        
        // Mettre à jour position
        _mobiles[i].p.x += _mobiles[i].v.x * dt;
        _mobiles[i].p.y += _mobiles[i].v.y * dt;
        
        // Collision avec les murs avec rebond
        if (_mobiles[i].p.x - _mobiles[i].r <= -1.0f) {
            _mobiles[i].v.x = -_mobiles[i].v.x * e;
            _mobiles[i].p.x = -1.0f + _mobiles[i].r;
        }
        if (_mobiles[i].p.x + _mobiles[i].r >= 1.0f) {
            _mobiles[i].v.x = -_mobiles[i].v.x * e;
            _mobiles[i].p.x = 1.0f - _mobiles[i].r;
        }
        if (_mobiles[i].p.y - _mobiles[i].r <= -1.0f) {
            _mobiles[i].v.y = -_mobiles[i].v.y * e;
            _mobiles[i].p.y = -1.0f + _mobiles[i].r;
        }
        if (_mobiles[i].p.y + _mobiles[i].r >= 1.0f) {
            _mobiles[i].v.y = -_mobiles[i].v.y * e;
            _mobiles[i].p.y = 1.0f - _mobiles[i].r;
        }
        
        // Amortissement global (facultatif)
        //_mobiles[i].v.x *= 0.95;//95f;
        //_mobiles[i].v.y *= 0.95;//95f;
        
        // Mettre à jour la couleur en fonction de la pression (visualisation)
        float pressure_ratio = (_mobiles[i].pressure / (GAS_CONSTANT * REST_DENSITY))*0.01f;
        pressure_ratio = fmaxf(0.0f, fminf(1.0f, pressure_ratio * 0.1f));
        
        _mobiles[i].color[0] = pressure_ratio;
        _mobiles[i].color[1] = 0.2f + 0.8f * (1.0f - pressure_ratio);
        _mobiles[i].color[2] = 1.0f - pressure_ratio;//TODO changer la couleur
		// Si la densité est trop élevée, réduire la vitesse
		//if (_mobiles[i].density > REST_DENSITY * 1.5f) {
		//	_mobiles[i].v.x *= 0.85f;  // Réduire davantage la vitesse dans les zones denses
		//	_mobiles[i].v.y *= 0.85f;
		//}
		/*
		// Calculer la vitesse actuelle
		float speed = sqrtf(_mobiles[i].v.x * _mobiles[i].v.x + _mobiles[i].v.y * _mobiles[i].v.y);
    
		// Si la vitesse dépasse le maximum, la réduire
		if (speed > max_speed) {
			float scale = max_speed / speed;
			_mobiles[i].v.x *= scale;
			_mobiles[i].v.y *= scale;
		}
		*/
    }
	for (int i = 0; i < _nb_mobiles; ++i) {
		for (int j = i + 1; j < _nb_mobiles; ++j) {
			float dx = _mobiles[j].p.x - _mobiles[i].p.x;
			float dy = _mobiles[j].p.y - _mobiles[i].p.y;
			float dist2 = dx * dx + dy * dy;
			
			// Somme des rayons
			float sumRadii = _mobiles[i].r + _mobiles[j].r;
			
			// Test de collision
			if (dist2 < sumRadii * sumRadii) {
				float dist = sqrtf(dist2);
				if (dist < 0.0001f) dist = 0.0001f; // éviter la division par zéro
				
				// Vecteur unitaire i -> j
				float nx = dx / dist;
				float ny = dy / dist;
				
				// Écarter les particules pour corriger le chevauchement
				float overlap = 0.5f * (sumRadii - dist);
				_mobiles[i].p.x -= nx * overlap;
				_mobiles[i].p.y -= ny * overlap;
				_mobiles[j].p.x += nx * overlap;
				_mobiles[j].p.y += ny * overlap;
				
				// Vitesse relative dans la direction de la collision
				float vx = _mobiles[j].v.x - _mobiles[i].v.x;
				float vy = _mobiles[j].v.y - _mobiles[i].v.y;
				float dot = vx * nx + vy * ny;
				
				// Si dot > 0, elles s'éloignent déjà ⇒ pas de correction supplémentaire
				if (dot > 0.0f) continue;
				
				// Coefficient de restitution (rebond)
				float restitution = 0.5f; // Ajustez selon l'effet rebond désiré
				
				// Impulsion
				float impulse = -(1.0f + restitution) * dot;
				// Masse = 1 pour les deux particules (dans votre code, MASS=1.0f)
				impulse *= 0.5f; // Répartition égale si les masses sont égales
				
				// Appliquer l’impulsion
				_mobiles[i].v.x -= nx * impulse;
				_mobiles[i].v.y -= ny * impulse;
				_mobiles[j].v.x += nx * impulse;
				_mobiles[j].v.y += ny * impulse;
			}
		}
	}
    rect_collide_all(_mobiles, _nb_mobiles, e);
}


// Fonction pour basculer entre les modes eau et espace
void toggle_simulation_mode() {
    SPACE_MODE = !SPACE_MODE; // Toggle the simulation mode
    mobile_quit();
    
    if (SPACE_MODE) {
        TIME_SCALE = 1000.0f;
        // Désactiver la gravité globale
        _g.x = _g.y = _g.z = 0.0f;
        // Initialiser en mode espace
        space_init(1023);
    } else {
        // Rétablir la gravité vers le bas
        _g.x = 0.0f; 
        _g.y = -9.81f;
        _g.z = 0.0f;
        // Réinitialiser en mode eau
        mobile_init(1023);
    }
}
// Fonction principale de simulation modifiée pour sélectionner le mode
void mobile_simu_with_mode_selection() {
    if (SPACE_MODE) {
        space_simulation();
    } else {
        mobile_simu();
    }
}

void mobile_draw(void){
	GLfloat *tmp = malloc(4 * _nb_mobiles * sizeof *tmp);
	assert(tmp);
	for (int i = 0; i < _nb_mobiles; ++i){
		tmp[4 * i + 0] = _mobiles[i].p.x;
		tmp[4 * i + 1] = _mobiles[i].p.y;
		tmp[4 * i + 2] = _mobiles[i].r;
	}
	glUniform4fv(glGetUniformLocation(_pId, "positions"), _nb_mobiles, tmp);
    /*
    if (SPACE_MODE) {
        //maj des couleurs en fonction de la vélocité ou autres propriétés
        for (int i = 0; i < _nb_mobiles; ++i) {
            //particule centrale (étoile) est plus brillante
            if (i == 0) {
                tmp[4 * i + 0] = 1.0f;       // Rouge plus vif
                tmp[4 * i + 1] = 0.9f;       // Jaune-orange
                tmp[4 * i + 2] = 0.3f;       // Un peu de bleu pour moins jaune
                tmp[4 * i + 3] = 1.0f;       // Opacité complète
            } else {
                // Copier les couleurs déjà calculées dans space_simulation
                for (int j = 0; j < 4; ++j) {
                    tmp[4 * i + j] = _mobiles[i].color[j];
                }
                
                // Ajouter un effet de brillance basé sur la vélocité
                float speed = sqrtf(_mobiles[i].v.x * _mobiles[i].v.x + 
                                   _mobiles[i].v.y * _mobiles[i].v.y +
                                   _mobiles[i].v.z * _mobiles[i].v.z);
                tmp[4 * i + 3] = 0.7f + 0.3f * speed; // Variation de l'opacité
            }
        }
    } else {// Mode eau normal : les couleurs de base
        for (int i = 0; i < _nb_mobiles; ++i) {
            for (int j = 0; j < 4; ++j) {
                tmp[4 * i + j] = _mobiles[i].color[j];
            }
        }
    }
    */
	for (int i = 0; i < _nb_mobiles; ++i)
		for (int j = 0; j < 4; ++j)
			tmp[4 * i + j] = _mobiles[i].color[j];
	glUniform4fv(glGetUniformLocation(_pId, "couleurs"), _nb_mobiles, tmp);
	glUniform1i(glGetUniformLocation(_pId, "nbe"), _nb_mobiles);
    glUniform1i(glGetUniformLocation(_pId, "space_mode"), SPACE_MODE);
	//TODO met le scatering et blur
	gl4dgDraw(_quad);
	free(tmp);
}

void mobile_quit(void) {
    if (_mobiles) {
        free(_mobiles);
        _mobiles = NULL;
        _nb_mobiles = 0;
    }
    
    if (_grid) {
        free(_grid);
        _grid = NULL;
    }
    
    if (_cell_indices) {
        free(_cell_indices);
        _cell_indices = NULL;
    }
    
    if (_particle_indices) {
        free(_particle_indices);
        _particle_indices = NULL;
    }
    rect_cleanup();
}
/*
void mobile_quit(void){
	if (_mobiles)
	{
		free(_mobiles);
		_mobiles = NULL;
		_nb_mobiles = 0;
	}
}
*/
