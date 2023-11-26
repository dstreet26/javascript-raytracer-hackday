const fs = require("fs");

const maxVal = 255;

const width = 255
const height = 255
const colors = [];

for (let i = 0; i < width; i++) {
      for (let j = 0; j < height; j++) {
        colors.push([i, j, 0]);
  }
}

fs.writeFileSync(
  "ppm-gradient.ppm",
  `P3\n${width} ${height}\n${maxVal}\n${
    colors.map((color) => color.join(" ")).join("\n") + "\n"
  }`
);
