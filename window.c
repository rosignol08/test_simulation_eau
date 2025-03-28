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
#include <math.h>

#define MAX_NEIGHBOURS 64
#define HASH_SIZE 1024
#define CELL_SIZE 0.1f  //taille des cellules pour la grille spatiale

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

//macro pour les expressions mathématiques complexes
#define POLY6 (315.0f / (64.0f * M_PI * powf(H, 9)))
#define SPIKY_GRAD (-45.0f / (M_PI * powf(H, 6)))
#define VISC_LAP (45.0f / (M_PI * powf(H, 6)))
#define SURFACE_TENSION 0.0728f

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
static float TIME_SCALE = 1.0f;  // 1.0 = vitesse normale, 0.5 = demi-vitesse


static void init(void);
static void draw(void);
static void quit(void);

static void mobile_init(int n);
static void mobile_simu(void);
static void mobile_draw(void);
static void mobile_quit(void);

/* on créé une variable pour stocker l'identifiant du programme GPU */
GLuint _pId = 0;

GLuint _quad = 0;

/* gravité */
static GLfloat _ig = 9.81f / 2.0f;
static vec3d_t _g = {0.0f, -98.1f, 0.0f}; // Modification ici: définir la gravité vers le bas à -9.81f
static const GLfloat e = 0.5f; //8.0f / 9.0f;

/* simulation d'eau de jsp qui */
// Ajouter ces paramètres SPH
static const float REST_DENSITY = 1000.0f;  // Densité au repos du fluide
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
        int cell_id = _mobiles[i].cell_id;
        
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
		    _mobiles[i].pressure = GAS_CONSTANT * 0.2f * (density_ratio - 1.0f);
		}
		//float density_ratio = _mobiles[i].density / REST_DENSITY;
        //_mobiles[i].pressure = GAS_CONSTANT * (powf(density_ratio, 7) - 1.0f);
        
        // Limiter les pressions extrêmes pour éviter l'instabilité
        if (_mobiles[i].pressure > GAS_CONSTANT * 10.0f)
            _mobiles[i].pressure = GAS_CONSTANT * 10.0f;
        if (_mobiles[i].pressure < -GAS_CONSTANT)
            _mobiles[i].pressure = -GAS_CONSTANT;
    }
    
    // Calculer les forces
    for (int i = 0; i < _nb_mobiles; i++) {
        _mobiles[i].force.x = 0.0f;
        _mobiles[i].force.y = 0.0f;
        _mobiles[i].force.z = 0.0f;
        
        int cell_id = _mobiles[i].cell_id;
        
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
                        
                        _mobiles[i].force.x += pressure_factor * dx / r;
                        _mobiles[i].force.y += pressure_factor * dy / r;
                        
                        // Force de viscosité
                        float visc_factor = VISCOSITY * MASS * 
                                         (_mobiles[j].v.x - _mobiles[i].v.x) * 
                                         //kernel_viscosity_laplacian(r) / _mobiles[j].density;
										 kernel_viscosity_improved(r, H);
                        
                        _mobiles[i].force.x += visc_factor;
                        
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
						        float repulsion_strength = 5000.0f * (min_distance * 1.1f - r) / min_distance;
						        _mobiles[i].force.x += repulsion_strength * dx / r;
						        _mobiles[i].force.y += repulsion_strength * dy / r;
						    } 
						    // Deuxième couche - répulsion plus douce
						    else {
						        float repulsion_strength = 500.0f * (min_distance * 1.5f - r) / min_distance;
						        _mobiles[i].force.x += repulsion_strength * dx / r;
						        _mobiles[i].force.y += repulsion_strength * dy / r;
						    }
						}
						//float min_distance = _mobiles[i].r + _mobiles[j].r;
						//if (r < min_distance) {
						//	float repulsion_strength = 1000.0f * (min_distance - r) / min_distance;
						//	_mobiles[i].force.x += repulsion_strength * dx / r;
						//	_mobiles[i].force.y += repulsion_strength * dy / r;
						//}
						
						// Limiter la densité
						if (_mobiles[i].density > REST_DENSITY * 2.0f) {
							_mobiles[i].v.x *= 0.5f;
							_mobiles[i].v.y *= 0.5f;
						}
						//if (r < min_distance * 1.5f) {
						//    float repulsion_strength = 200.0f * (min_distance * 1.5f - r) / min_distance;
						//    _mobiles[i].force.x -= repulsion_strength * dx / r;
						//    _mobiles[i].force.y -= repulsion_strength * dy / r;
						//}
                    }
                }
            }
        }
        
        // Ajouter la gravité
        //_mobiles[i].force.x += _g.x * _mobiles[i].density * 1.0f;
        //_mobiles[i].force.y += _g.y * _mobiles[i].density * 1.0f;
		//limiter leur magnitude
        float force_magnitude = sqrtf(_mobiles[i].force.x * _mobiles[i].force.x + 
                                     _mobiles[i].force.y * _mobiles[i].force.y);
        const float max_force = 5000.0f;
        
        if (force_magnitude > max_force) {
            float scale = max_force / force_magnitude;
            _mobiles[i].force.x *= scale;
            _mobiles[i].force.y *= scale;
        }
		//_mobiles[i].force.x += _g.x * _mobiles[i].density * 0.20f;
		//_mobiles[i].force.y += _g.y * _mobiles[i].density * 0.20f;
		_mobiles[i].force.x += _g.x * MASS * 100.0f;  // Indépendant de la densité
		_mobiles[i].force.y += _g.y * MASS * 100.0f;  // Indépendant de la densité

    }
}
/*fin fonction de ouf*/

/* initialise des paramètres GL et GL4D */
void init(void){
	_quad = gl4dgGenQuadf();
	/* activer la synchronisation verticale */
	SDL_GL_SetSwapInterval(1);
	/* set la couleur d'effacement OpenGL */
	glClearColor(0.0f, 0.0f, 0.5f, 1.0f);
	/* créer un programme GPU pour OpenGL (en GL4D) */
	_pId = gl4duCreateProgram("<vs>shaders/identity.vs", "<fs>shaders/calculs.fs", NULL);

	mobile_init(1024);
}

void draw(void){
	/* effacer le buffer de couleur (image) et le buffer de profondeur d'OpenGL */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* utiliser le programme GPU "_pId" */
	glUseProgram(_pId);
	/* binder (mettre au premier plan, "en courante" ou "en active") la
	   matrice view */

	mobile_draw();

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
    
    // Si on n'a pas assez de particules, remplir le reste
    for (int i = index; i < n; i++) {
        _mobiles[i].p.x = gl4dmSURand() * 0.5f - 0.5f;
        _mobiles[i].p.y = gl4dmSURand() * 0.5f + 0.5f;
        _mobiles[i].p.z = 0.0f;
        _mobiles[i].v.x = 0.0f;
        _mobiles[i].v.y = 0.0f;
        _mobiles[i].v.z = 0.0f;
        _mobiles[i].r = 0.02f;
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
	/*
	assert(_mobiles == NULL);
	_nb_mobiles = n;
	_mobiles = malloc(_nb_mobiles * sizeof *_mobiles);
	assert(_mobiles);
	for (int i = 0; i < _nb_mobiles; ++i){
		_mobiles[i].p.x = 0.89f * gl4dmSURand();
		_mobiles[i].p.y = 0.89f * gl4dmSURand() + 0.1f; // spawn les balles plus haut
		_mobiles[i].p.z = 0.0f;							// 0.89f * gl4dmSURand(); /* Ajout de la coordonnée z */
//		_mobiles[i].v.x = gl4dmSURand() * 0.5f;			// plus rapide au début
//		_mobiles[i].v.y = gl4dmSURand() * 0.5f;							// gl4dmSURand() * 0.5f; //plus rapide au début
//		_mobiles[i].v.z = 0.0f; 			/* Vitesse initiale en z */
//		_mobiles[i].r = 0.02f; //taille uniforme //0.01f + 0.1f * gl4dmURand();
//		_mobiles[i].color[0] = gl4dmURand();
//		_mobiles[i].color[1] = gl4dmURand();
//		_mobiles[i].color[2] = gl4dmURand();
//		_mobiles[i].color[3] = 1.0f;
//	}
	
}

void mobile_simu(void){
	static double t0 = 0;
	double t = gl4dGetElapsedTime() / 1000.0, dt = (t - t0) * TIME_SCALE;
	//double t = gl4dGetElapsedTime() / 1000.0, dt = t - t0;
	t0 = t;

	if (dt > 0.03f) dt = 0.03f; // Limiter le pas de temps à 30 ms
	//const float max_speed = 0.5f; // Vitesse maximale autorisée

	// Calculer les forces SPH
	compute_sph_forces();
	for (int i = 0; i < _nb_mobiles; ++i) {
        // Intégration explicite d'Euler
        float accel_x = _mobiles[i].force.x / _mobiles[i].density;
        float accel_y = _mobiles[i].force.y / _mobiles[i].density;
        
        _mobiles[i].v.x += accel_x * dt;
        _mobiles[i].v.y += accel_y * dt;
        
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
        float pressure_ratio = _mobiles[i].pressure / (GAS_CONSTANT * REST_DENSITY);
        pressure_ratio = fmaxf(0.0f, fminf(1.0f, pressure_ratio * 0.1f));
        
        _mobiles[i].color[0] = pressure_ratio;
        _mobiles[i].color[1] = 0.2f + 0.8f * (1.0f - pressure_ratio);
        _mobiles[i].color[2] = 1.0f - pressure_ratio;
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
}
	/*
	const float seuil_collision = 0.05f; // seuil de vitesse pour la collision

	for (int i = 0; i < _nb_mobiles; ++i){
		int collision_sol = 0, collision = 0;
		
		_mobiles[i].p.x += _mobiles[i].v.x * dt;
		_mobiles[i].p.y += _mobiles[i].v.y * dt;
		//_mobiles[i].p.z += _mobiles[i].v.z * dt; //pour que ça soit en 2D

		//collisions avec les murs en x
		if (_mobiles[i].p.x - _mobiles[i].r <= -1.0f){
			if (_mobiles[i].v.x < 0.0f)
				_mobiles[i].v.x = -_mobiles[i].v.x;
				_mobiles[i].p.x = -1.0f + _mobiles[i].r; // Repositionner
			collision = 1;
		}
		if (_mobiles[i].p.x + _mobiles[i].r >= 1.0f){
			if (_mobiles[i].v.x > 0.0f)
				_mobiles[i].v.x = -_mobiles[i].v.x;
				_mobiles[i].p.x = 1.0f - _mobiles[i].r; // Repositionner
			collision = 1;
		}
		
		////collision avec les murs en z
		//if (_mobiles[i].p.z - _mobiles[i].r <= -1.0f){
		//	if (_mobiles[i].v.z < 0.0f)
		//		_mobiles[i].v.z = -_mobiles[i].v.z;
		//		_mobiles[i].p.z = -1.0f + _mobiles[i].r; // Repositionner
		//	collision = 1;
		//}
		//if (_mobiles[i].p.z + _mobiles[i].r >= 1.0f){
		//	if (_mobiles[i].v.z > 0.0f)
		//		_mobiles[i].v.z = -_mobiles[i].v.z;
		//		_mobiles[i].p.z = 1.0f - _mobiles[i].r; // Repositionner
		//	collision = 1;
		//}
		
		//collision avec le sol
		if (_mobiles[i].p.y - _mobiles[i].r <= -1.0f){
			if (_mobiles[i].v.y < 0.0f)
				_mobiles[i].v.y = -_mobiles[i].v.y;
				_mobiles[i].p.y = -1.0f + _mobiles[i].r; // Repositionner
			collision = 1;
			collision_sol = 1;

			// si la balle touche le sol et sa vitesse est inférieure au seuil
			// on la stoppe
			if (fabs(_mobiles[i].v.y) < seuil_collision){
				_mobiles[i].v.y = 0.0f;
				//_mobiles[i].p.y = -1.0f + _mobiles[i].r; // repositionne la balle au dessus du sol TODO changer ça pour les colision entre balles
			}
		}
		//collision avec le plafond
		if (_mobiles[i].p.y + _mobiles[i].r >= 1.0f){
			if (_mobiles[i].v.y > 0.0f)
				_mobiles[i].v.y = -_mobiles[i].v.y;
				_mobiles[i].p.y = 1.0f - _mobiles[i].r; // Repositionner
			collision = 1;
		}

		//application d'amortissement si collision
		if (collision != 0){
			_mobiles[i].v.x *= e;
			_mobiles[i].v.y *= e;
			_mobiles[i].v.z *= e;
			
			//si la vitesse est inférieure au seuil, on la met à 0
			if(fabs(_mobiles[i].v.x) < seuil_collision)
    		  _mobiles[i].v.x = 0.0f;
    		if(fabs(_mobiles[i].v.y) < seuil_collision)
    		  _mobiles[i].v.y = 0.0f;
    		//if(fabs(_mobiles[i].v.z) < seuil_collision)
    		//  _mobiles[i].v.z = 0.0f;
    	}
		if (collision_sol == 0){
			_mobiles[i].v.x += _g.x * dt;
			_mobiles[i].v.y += _g.y * dt;
			//_mobiles[i].v.z += _g.z * dt;
		}
		//// Si la balle est au sol et n'a plus de vitesse horizontale significative,
    	//// elle ne devrait plus bouger horizontalement
    	//else if(_mobiles[i].v.y == 0.0f) {
		//	_mobiles[i].v.x = 0.0f;
		//	_mobiles[i].v.z = 0.0f;
		//}
		
		// Si la balle est au sol et n'a plus de vitesse verticale significative
		if (collision_sol && _mobiles[i].v.y == 0.0f) {
			// Réduire progressivement la vitesse horizontale (friction)
			_mobiles[i].v.x *= 0.98f;
			//_mobiles[i].v.z *= 0.98f;

			// Arrêter complètement si très lent
			if (fabs(_mobiles[i].v.x) < seuil_collision)
				_mobiles[i].v.x = 0.0f;
			//if (fabs(_mobiles[i].v.z) < seuil_collision)
			//	_mobiles[i].v.z = 0.0f;
		}
	}
	//gestion des collisions entre balles
    for (int i = 0; i < _nb_mobiles; ++i) {
        for (int j = i + 1; j < _nb_mobiles; ++j) {
            // Calculer la distance entre les centres des deux balles
            float dx = _mobiles[j].p.x - _mobiles[i].p.x;
            float dy = _mobiles[j].p.y - _mobiles[i].p.y;
            //float dz = _mobiles[j].p.z - _mobiles[i].p.z;
            
            // Distance au carré (évite de calculer sqrt pour l'optimisation)
            float distanceSquared = dx * dx + dy * dy; //+ dz * dz;
            
            // Somme des rayons
            float sumRadii = _mobiles[i].r + _mobiles[j].r;
            
            // Test de collision - collision si la distance est inférieure à la somme des rayons
            if (distanceSquared < sumRadii * sumRadii) {
                // Calcul de la distance réelle
                float distance = sqrt(distanceSquared);
                
                // Éviter division par zéro
                if (distance == 0.0f) distance = 0.0001f;
                
                // Vecteur unitaire de la direction de collision (de i vers j)
                float nx = dx / distance;
                float ny = dy / distance;
                //float nz = dz / distance;
                
                // Calculer le chevauchement pour repositionner les balles
                float overlap = (sumRadii - distance) / 2.0f;
                
                // Repositionner les balles
                _mobiles[i].p.x -= nx * overlap;
                _mobiles[i].p.y -= ny * overlap;
                //_mobiles[i].p.z -= nz * overlap;
                
                _mobiles[j].p.x += nx * overlap;
                _mobiles[j].p.y += ny * overlap;
                //_mobiles[j].p.z += nz * overlap;
                
                // Calcul de la vitesse relative dans la direction de collision
                float vx = _mobiles[j].v.x - _mobiles[i].v.x;
                float vy = _mobiles[j].v.y - _mobiles[i].v.y;
                //float vz = _mobiles[j].v.z - _mobiles[i].v.z;
                
                // Projection de la vitesse relative sur la normale de collision
                float dotProduct = nx * vx + ny * vy; //+ nz * vz;
                
                // Si les balles s'éloignent déjà, pas besoin de les rebondir
                if (dotProduct >= 0.0f) continue;
                
                // Coefficient de restitution (rebond)
                float restitution = e;
                
                // Calculer l'impulsion
                float impulse = -(1.0f + restitution) * dotProduct;
                
                // Masse égale pour toutes les balles (on peut ajuster ici si nécessaire)
                impulse /= 2.0f;
                
                // Appliquer l'impulsion aux vitesses
                _mobiles[i].v.x -= nx * impulse;
                _mobiles[i].v.y -= ny * impulse;
                //_mobiles[i].v.z -= nz * impulse;
                
                _mobiles[j].v.x += nx * impulse;
                _mobiles[j].v.y += ny * impulse;
                //_mobiles[j].v.z += nz * impulse;
            }
        }
	}
}
*/
void mobile_draw(void){
	GLfloat *tmp = malloc(4 * _nb_mobiles * sizeof *tmp);
	assert(tmp);
	for (int i = 0; i < _nb_mobiles; ++i)
	{
		tmp[4 * i + 0] = _mobiles[i].p.x;
		tmp[4 * i + 1] = _mobiles[i].p.y;
		tmp[4 * i + 2] = _mobiles[i].r;
	}
	glUniform4fv(glGetUniformLocation(_pId, "positions"), _nb_mobiles, tmp);
	for (int i = 0; i < _nb_mobiles; ++i)
		for (int j = 0; j < 4; ++j)
			tmp[4 * i + j] = _mobiles[i].color[j];
	glUniform4fv(glGetUniformLocation(_pId, "couleurs"), _nb_mobiles, tmp);
	glUniform1i(glGetUniformLocation(_pId, "nbe"), _nb_mobiles);
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