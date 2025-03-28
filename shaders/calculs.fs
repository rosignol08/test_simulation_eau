#version 330
out vec4 fragColor;
in vec2 fcoord;
uniform vec4 positions[1024];
uniform vec4 couleurs[1024];
uniform int nbe;

uniform vec4 rectangles[1024]; // Rectangles : (x, y, largeur, hauteur)
uniform vec4 rect_color;       // Couleur des rectangles
uniform int nb_rects;          // Nombre de rectangles

void main() {
  //fragColor = vec4(gl_FragCoord.xy / 600.0, 0.0, 1.0);
  //si le fragment est dans un rectangle
  for (int i = 0; i < nb_rects; i++) {
    vec4 rect = rectangles[i];
    if (fcoord.x >= rect.x && fcoord.x <= rect.x + rect.z &&
      fcoord.y >= rect.y && fcoord.y <= rect.y + rect.w) {
      fragColor = rect_color;
      return;
    }
  }
  //sinon 
  for (int i = 0; i < nbe; ++i) {
    if(distance(positions[i].xy, fcoord) < positions[i].z /* le rayon du mobile i */) {
      fragColor = couleurs[i];
      return;
    }      
  }
  fragColor = vec4(0.0, 0.0, 1.0, 1.0);
}