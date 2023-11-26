const fs = require("fs");

const maxVal = 255;

const colors = [
  [255, 0, 0],
  [0, 255, 0],
  [0, 0, 255],
  [255, 255, 0],
  [255, 0, 255],
  [0, 255, 255],
];

const width = 2;
const height = 3;

for (let i = 0; i < width; i++) {
  for (let j = 0; j < height; j++) {
    colors.push([i, j, 0]);
  }
}

fs.writeFileSync(
  "sixpixels.ppm",
  `P3\n${width} ${height}\n${maxVal}\n${
    colors.map((color) => color.join(" ")).join("\n") + "\n"
  }`
);
