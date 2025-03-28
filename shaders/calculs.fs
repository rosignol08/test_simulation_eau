#version 330
out vec4 fragColor;
in vec2 fcoord;
uniform vec4 positions[1024];
uniform vec4 couleurs[1024];
uniform int nbe;

uniform vec4 rectangles[10]; // Rectangles : (x, y, largeur, hauteur)
uniform vec4 rect_color;       // Couleur des rectangles
uniform int nb_rects;          // Nombre de rectangles
uniform float rect_angles[10]; // Angles des rectangles

void main() {
  //fragColor = vec4(gl_FragCoord.xy / 600.0, 0.0, 1.0);
  //si le fragment est dans un rectangle
  for (int i = 0; i < nb_rects; i++) {
    vec4 rect = rectangles[i];
    float angle = rect_angles[i];
    vec2 center = rect.xy;// + vec2(rect.z / 2.0, rect.w / 2.0); // centre du rectangle
    // Translate the fragment coordinate to the rectangle's local space
    vec2 localCoord = fcoord - center;

    // Apply rotation
    float cosAngle = cos(angle);
    float sinAngle = sin(angle);
    vec2 rotatedCoord = vec2(
      cosAngle * localCoord.x + sinAngle * localCoord.y,
      -sinAngle * localCoord.x + cosAngle * localCoord.y
    );

    // Check if the rotated coordinate is within the rectangle's bounds
    if (rotatedCoord.x >= 0.0 && rotatedCoord.x <= rect.z &&
      rotatedCoord.y >= 0.0 && rotatedCoord.y <= rect.w) {
      
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