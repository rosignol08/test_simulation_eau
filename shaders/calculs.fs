#version 330
out vec4 fragColor;
in vec2 fcoord;
uniform vec4 positions[256];
uniform vec4 couleurs[256];
uniform int nbe;

void main() {
  //fragColor = vec4(gl_FragCoord.xy / 600.0, 0.0, 1.0);
  for(int i = 0; i < nbe; ++i) {
    if(distance(positions[i].xy, fcoord) < positions[i].z /* le rayon du mobile i */) {
      fragColor = couleurs[i];
      return;
    }      
  }
  fragColor = vec4(0.0, 0.0, 1.0, 1.0);
}
