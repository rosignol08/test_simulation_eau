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

typedef struct vec3d_t vec3d_t;
typedef struct mobile_t mobile_t;

struct vec3d_t
{
	GLfloat x, y, z;
};

struct mobile_t
{
	vec3d_t p, v;
	GLfloat r;
	GLfloat color[4];
};

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
static vec3d_t _g = {0.0f, -1.81f, 0.0f}; // Modification ici: définir la gravité vers le bas à -9.81f
static const GLfloat e = 8.0f / 9.0f;

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

/* initialise des paramètres GL et GL4D */
void init(void){
	_quad = gl4dgGenQuadf();
	/* activer la synchronisation verticale */
	SDL_GL_SetSwapInterval(1);
	/* set la couleur d'effacement OpenGL */
	glClearColor(0.0f, 0.0f, 0.5f, 1.0f);
	/* créer un programme GPU pour OpenGL (en GL4D) */
	_pId = gl4duCreateProgram("<vs>shaders/identity.vs", "<fs>shaders/calculs.fs", NULL);

	mobile_init(230);
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
	for (int i = 0; i < _nb_mobiles; ++i){
		_mobiles[i].p.x = 0.89f * gl4dmSURand();
		_mobiles[i].p.y = 0.89f * gl4dmSURand() + 0.1f; // spawn les balles plus haut
		_mobiles[i].p.z = 0.0f;							// 0.89f * gl4dmSURand(); /* Ajout de la coordonnée z */
		_mobiles[i].v.x = gl4dmSURand() * 0.5f;			// plus rapide au début
		_mobiles[i].v.y = gl4dmSURand() * 0.5f;							// gl4dmSURand() * 0.5f; //plus rapide au début
		_mobiles[i].v.z = 0.0f; 			/* Vitesse initiale en z */
		_mobiles[i].r = 0.02f; //taille uniforme //0.01f + 0.1f * gl4dmURand();
		_mobiles[i].color[0] = gl4dmURand();
		_mobiles[i].color[1] = gl4dmURand();
		_mobiles[i].color[2] = gl4dmURand();
		_mobiles[i].color[3] = 1.0f;
	}
}

void mobile_simu(void){
	static double t0 = 0;
	double t = gl4dGetElapsedTime() / 1000.0, dt = t - t0;
	t0 = t;
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
		/*
		//collision avec les murs en z
		if (_mobiles[i].p.z - _mobiles[i].r <= -1.0f){
			if (_mobiles[i].v.z < 0.0f)
				_mobiles[i].v.z = -_mobiles[i].v.z;
				_mobiles[i].p.z = -1.0f + _mobiles[i].r; // Repositionner
			collision = 1;
		}
		if (_mobiles[i].p.z + _mobiles[i].r >= 1.0f){
			if (_mobiles[i].v.z > 0.0f)
				_mobiles[i].v.z = -_mobiles[i].v.z;
				_mobiles[i].p.z = 1.0f - _mobiles[i].r; // Repositionner
			collision = 1;
		}
		*/
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
		/*
		// Si la balle est au sol et n'a plus de vitesse horizontale significative,
    	// elle ne devrait plus bouger horizontalement
    	else if(_mobiles[i].v.y == 0.0f) {
			_mobiles[i].v.x = 0.0f;
			_mobiles[i].v.z = 0.0f;
		}
		*/
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
}

void mobile_quit(void){
	if (_mobiles)
	{
		free(_mobiles);
		_mobiles = NULL;
		_nb_mobiles = 0;
	}
}
