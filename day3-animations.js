// node day3-animations.js
// ffmpeg  -i 'animations/anim%d.png' -c:v libx264 -pix_fmt yuv420p out.mp4
// ffplay out.mp4

const Jimp = require("jimp");
const imageWidth = 100;
const animationFrames = 100;
const aspectRatio = 16 / 9;
const imageHeight = Math.floor(imageWidth / aspectRatio);
const fov = 50;
const focusDist = 20;
const samplesPerPixel = 10;
const maxDepth = 10;
const viewportHeight = 2 * Math.tan((fov * (Math.PI / 180)) / 2) * focusDist;
const viewportWidth = viewportHeight * (imageWidth / imageHeight);
const cameraCenter = [13, 2, 3];
const lookAt = [0, 0, 0];
const vUp = [0, 1, 0];
const w = unitVector(subtractVectors(cameraCenter, lookAt));
const u = unitVector(crossVectors(vUp, w));
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
const pixel00Location = addVectors(
  viewportUpperLeft,
  multiplyVector(addVectors(pixelDeltaU, pixelDeltaV), 0.5)
);
const maxVal = 255;

//range foreach because we need to await on each image creation
(async () => {
  const range = (start, stop, step) =>
    Array.from(
      { length: (stop - start) / step + 1 },
      (_, i) => start + i * step
    );
  for (const t of range(0, animationFrames, 1)) {
    await createImage(t);
  }
})();

function createWorld(t) {
  const metalSphere1 = {
    type: "sphere",
    center: [4, 1, 4],
    radius: 1,
    material: {
      type: "metal",
      albedo: [0.8, 0.8, 0],
      fuzz: 0,
    },
  };
  const metalSphere2 = {
    type: "sphere",
    center: [4, 1, 6],
    radius: 1,
    material: {
      type: "metal",
      albedo: [0.8, 0, 0.8],
      fuzz: 0.0,
    },
  };
  const bigSphere = {
    type: "sphere",
    center: [0, -100, -1],
    radius: 100,
    material: {
      type: "lambertian",
      albedo: [0.5, 0.5, 0.5],
    },
  };
  const sphere1 = {
    type: "sphere",
    // center: [6, 0, (-t*0.2)+7],
    center: [4, 0, (t * -0.2) + 10],
    radius: 2,
    material: {
      type: "glass",
      albedo: [1, 1, 1],
      ir: 1.2,
    },
  };
  const sphere2 = {
    type: "sphere",
    center: [0, 0, -2],
    radius: 2,
    material: {
      type: "lambertian",
      albedo: [0, 0, 1],
    },
  };
  const world = [sphere1, sphere2, bigSphere, metalSphere1, metalSphere2];
  return world;
}

async function createImage(t, callback) {
  // console.log("createImage", t);
  var world = createWorld(t);
  // console.log("world", world[0].center);
  // var s = fs.createWriteStream("image15.ppm");
  // s.write(`P3\n${imageWidth} ${imageHeight}\n${maxVal}\n`);

  var colors = [];
  for (let j = 0; j < imageHeight; j++) {
    // console.log(`${t}/${j}/${imageHeight}`)
    // process.stdout.clearLine(0);
    process.stdout.cursorTo(0);
    process.stdout.write(`${j}/${imageHeight} - ${t}/${animationFrames}`);
    for (let i = 0; i < imageWidth; i++) {
      let color = [0, 0, 0];
      const r = getRay(i, j);
      for (let k = 0; k < samplesPerPixel; k++) {
        color = addVectors(color, rayColor(r, maxDepth, world));
      }
      colors.push(writeColor(color, samplesPerPixel));
    }
  }

  return new Promise((resolve, reject) => {
    new Jimp(imageWidth, imageHeight, async function (err, image) {
      var data = new Uint8Array(imageWidth * imageHeight * 4);
      for (let i = 0; i < imageWidth * imageHeight; i++) {
        data[i * 4] = colors[i][0];
        data[i * 4 + 1] = colors[i][1];
        data[i * 4 + 2] = colors[i][2];
        data[i * 4 + 3] = maxVal;
      }
      this.bitmap.data = data;
      // console.log("animations/anim_0" + t + ".png");
      // await this.writeAsync("animations/anim_0" + t + ".png");
      await this.writeAsync("animations/anim" + t + ".png");
      resolve();
    });
  });
}

// s.end();
function getRay(i, j) {
  const pixelCenter = addVectors(
    pixel00Location,
    addVectors(multiplyVector(pixelDeltaU, i), multiplyVector(pixelDeltaV, j))
  );
  const rayDirection = subtractVectors(pixelCenter, cameraCenter);
  return [cameraCenter, rayDirection];
}
function rayColor(ray, depth, world) {
  if (depth <= 0) {
    return [0, 0, 0];
  }
  let min1 = 0.001;
  let max1 = 100000;
  let rec = null;
  for (let k = 0; k < world.length; k++) {
    const element = world[k];
    if (element.type === "sphere") {
      const hit = hitSphere(element.center, element.radius, ray, min1, max1);
      if (hit) {
        hit_anything = true;
        max1 = hit.t;
        rec = {
          ...hit,
          material: element.material,
        };
      }
    }
  }
  if (rec) {
    const scattered = scatter(ray, rec);
    if (scattered) {
      return multiplyVectors(
        rayColor(scattered.scatteredRay, depth - 1, world),
        scattered.attenuation
      );
    }
    return [0, 0, 0];
  }
  const a = 0.5 * (unitVector(ray[1])[1] + 1.0);
  return addVectors(
    multiplyVector([1.0, 1.0, 1.0], 1.0 - a),
    multiplyVector([0.5, 0.7, 1.0], a)
  );
}
function scatter(ray, rec) {
  if (rec.material.type === "lambertian") {
    return {
      scatteredRay: [rec.p, addVectors(rec.normal, randomUnitVector())],
      attenuation: rec.material.albedo,
    };
  } else if (rec.material.type === "metal") {
    const reflected = reflect(unitVector(ray[1]), rec.normal);
    const scattered = [
      rec.p,
      addVectors(
        reflected,
        multiplyVector(randomUnitVector(), rec.material.fuzz)
      ),
    ];
    if (dotVectors(scattered[1], rec.normal) > 0) {
      return {
        scatteredRay: scattered,
        attenuation: rec.material.albedo,
      };
    }
  } else if (rec.material.type === "glass") {
    const refractionRatio = rec.frontFace
      ? 1 / rec.material.ir
      : rec.material.ir;
    const unitDirection = unitVector(ray[1]);
    const cosTheta = Math.min(
      dotVectors(multiplyVector(unitDirection, -1), rec.normal),
      1
    );
    const sinTheta = Math.sqrt(1.0 - cosTheta * cosTheta);
    const cannotRefract = refractionRatio * sinTheta > 1.0;
    let direction = null;
    if (cannotRefract || reflectance(cosTheta, refractionRatio) > 0.5) {
      direction = reflect(unitDirection, rec.normal);
    } else {
      direction = refract(unitDirection, rec.normal, refractionRatio);
    }
    return {
      scatteredRay: [rec.p, direction],
      attenuation: rec.material.albedo,
    };
  }
  return null;
}
function hitSphere(center, radius, ray, min, max) {
  //todo: benchmark difference between dot and length_squared
  const oc = subtractVectors(ray[0], center);
  const a = lengthSquared(ray[1]);
  const halfB = dotVectors(oc, ray[1]);
  const c = lengthSquared(oc) - radius * radius;
  const discriminant = halfB * halfB - a * c;
  if (discriminant < 0) {
    return null;
  } else {
    const sqrtd = Math.sqrt(discriminant);
    let root = (-halfB - sqrtd) / a;
    if (!surrounds(root, min, max)) {
      root = (-halfB + sqrtd) / a;
      if (!surrounds(root, min, max)) {
        return null;
      }
    }
    const p = addVectors(ray[0], multiplyVector(ray[1], root));
    const outwardNormal = multiplyVector(
      subtractVectors(p, center),
      1 / radius
    );
    const frontFace = dotVectors(ray[1], outwardNormal) < 0;
    const normal = frontFace
      ? outwardNormal
      : multiplyVector(outwardNormal, -1);
    return {
      t: root,
      p,
      frontFace: frontFace,
      normal,
    };
  }
}
function subtractVectors(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}
function addVectors(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}
function multiplyVector(a, b) {
  return [a[0] * b, a[1] * b, a[2] * b];
}
function multiplyVectors(a, b) {
  return [a[0] * b[0], a[1] * b[1], a[2] * b[2]];
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
function lengthSquared(v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
function surrounds(x, min, max) {
  return min <= x && x <= max;
}
function random_in_unit_sphere() {
  while (true) {
    const p = randomVector3(-1, 1);
    if (lengthSquared(p) < 1) {
      return p;
    }
  }
}
function randomVector3(min, max) {
  return [randomFloat(min, max), randomFloat(min, max), randomFloat(min, max)];
}
function randomFloat(min, max) {
  return Math.random() * (max - min) + min;
}
function randomUnitVector() {
  return unitVector(random_in_unit_sphere());
}
function reflect(v, n) {
  return subtractVectors(v, multiplyVector(n, 2 * dotVectors(v, n)));
}
function refract(uv, n, etai_over_etat) {
  const cos_theta = Math.min(dotVectors(multiplyVector(uv, -1), n), 1.0);
  const r_out_perp = multiplyVector(
    addVectors(uv, multiplyVector(n, cos_theta)),
    etai_over_etat
  );
  const r_out_parallel = multiplyVector(
    n,
    -Math.sqrt(Math.abs(1.0 - lengthSquared(r_out_perp)))
  );
  return addVectors(r_out_perp, r_out_parallel);
}
function reflectance(cosine, ref_idx) {
  let r0 = (1 - ref_idx) / (1 + ref_idx);
  r0 = r0 * r0;
  return r0 + (1 - r0) * Math.pow(1 - cosine, 5);
}
function writeColor(color, samples_per_pixel) {
  const scale = 1.0 / samples_per_pixel;
  // const scaledColor = multiplyVector(color, scale);
  const scaledColor = multiplyVector(color, scale).map((x) => Math.sqrt(x));
  // const scaledColor2 = multiplyVector(scaledColor, 1);
  const colorScaled = multiplyVector(scaledColor, maxVal).map((x) =>
    Math.round(x).toFixed(0)
  );
  return colorScaled;
  //   s.write(`${colorScaled.join(" ")}\n`);
}
