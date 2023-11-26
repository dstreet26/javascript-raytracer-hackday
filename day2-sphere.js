const fs = require("fs");
const imageWidth = 200;
const aspectRatio = 16 / 9;
const imageHeight = Math.floor(imageWidth / aspectRatio);
const fov = 40;
const focusDist = 10;
const theta = degrees_to_radians(fov);
const h = Math.tan(theta / 2);
const viewportHeight = 2 * h * focusDist;
const viewportWidth = viewportHeight * (imageWidth / imageHeight);
const cameraCenter = [5, 5, 5];
const lookat = [0, 0, 0];
const vup = [0, 1, 0];
const w = unitVector(subtractVectors(cameraCenter, lookat));
const u = unitVector(crossVectors(vup, w));
const v = crossVectors(w, u);
const viewportU = multiplyVector(u, viewportWidth);
const viewportV = multiplyVector(v, -viewportHeight);
const pixelDeltaU = multiplyVector(viewportU, 1 / imageWidth);
const pixelDeltaV = multiplyVector(viewportV, 1 / imageHeight);

const viewportUpperLeft = subtractVectors(
  subtractVectors(
    subtractVectors(cameraCenter, multiplyVector(w, focusDist)),
    multiplyVector(viewportU, 0.5)
  ),
  multiplyVector(viewportV, 0.5)
);
// note: first iteration, broken and i can't explain why
// const viewport_upper_left_broken = subtractVectors(
//   camera_center,
//   subtractVectors(
//     multiplyVector(w, focus_dist),
//     subtractVectors(
//       multiplyVector(viewport_u, 0.5),
//       multiplyVector(viewport_v, 0.5)
//     )
//   )
// );

const pixel00Location = addVectors(
  viewportUpperLeft,
  multiplyVector(addVectors(pixelDeltaU, pixelDeltaV), 0.5)
);

const maxVal = 255;
var s = fs.createWriteStream("sphere.ppm");
s.write(`P3\n${imageWidth} ${imageHeight}\n${maxVal}\n`);

for (let j = 0; j < imageHeight; j++) {
  for (let i = 0; i < imageWidth; i++) {
    const pixelCenter = addVectors(
      pixel00Location,
      addVectors(multiplyVector(pixelDeltaU, i), multiplyVector(pixelDeltaV, j))
    );
    const rayDirection = subtractVectors(pixelCenter, cameraCenter);
    const unitDirection = unitVector(rayDirection);

    let color = multiplyVector(addVector(unitDirection, 1.0), 1.5);
    const t = hitSphere([0, 0, 0], 1, [cameraCenter, rayDirection]);
    if (t > 0.0) {
      const N = unitVector(
        subtractVectors(
          addVectors(cameraCenter, multiplyVector(rayDirection, t)),
          [0, 0, -1]
        )
      );
      color = multiplyVector([N[0] + 1, N[1] + 1, N[2] + 1], 0.5);
    }
    const colorScaled = multiplyVector(color, maxVal).map((x) =>
      Math.round(x).toFixed(0)
    );
    s.write(`${colorScaled.join(" ")}\n`);
  }
}

s.end();

function hitSphere(center, radius, ray) {
  //todo: benchmark difference between dot and length_squared
  const oc = subtractVectors(ray[0], center);
  // const a = dotVectors(ray[1], ray[1]);
  const a = length_squared(ray[1]);
  // const b = 2.0 * dotVectors(oc, ray[1]);
  const half_b = dotVectors(oc, ray[1]);
  // const c = dotVectors(oc, oc) - radius * radius;
  const c = length_squared(oc) - radius * radius;
  // const discriminant = b * b - 4 * a * c;
  const discriminant = half_b * half_b - a * c;

  if (discriminant < 0) {
    return -1.0;
  } else {
    //   return (-b - Math.sqrt(discriminant)) / (2.0 * a);
    return (-half_b - Math.sqrt(discriminant)) / a;
  }
}

function subtractVectors(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}
function addVectors(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}
function addVector(a, b) {
  return [a[0] + b, a[1] + b, a[2] + b];
}
function multiplyVector(a, b) {
  return [a[0] * b, a[1] * b, a[2] * b];
}
function unitVector(vector) {
  return multiplyVector(vector, 1 / vectorLength(vector));
}
function vectorLength(vector) {
  //todo: benchmark difference between hypot and manual
  return Math.hypot(vector[0], vector[1], vector[2]);
}
function dotVectors(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
function crossVectors(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    -(a[0] * b[2] - a[2] * b[0]),
    a[0] * b[1] - a[1] * b[0],
  ];
}
function degrees_to_radians(degrees) {
  var pi = Math.PI;
  return degrees * (pi / 180);
}
function length_squared(v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
