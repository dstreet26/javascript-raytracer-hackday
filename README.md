# Ray Tracer Hackday 

Mostly following the book [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html) by Peter Shirley.

However I did not want to use Classes or Operator/Prototype Overloading. It wasn't always clear if a vector was being multiplied by a scalar or another vector. I also wanted to avoid using a library for vectors and matrices but also looked into how some of them were implemented.

This video also follows the book but in Rust and I agree with a lot of the points made about the book:
https://youtu.be/6D8WVYm1YwY?si=ywrMFzKMB35x2wBh

Some other issues that I had with the book is that there is no light source or shadows in the first section, which I thought were fundamental to ray tracers.

Also, there was a part in the book about using an "interval", but the reason for needing it didn't seem to be explained well. The values for it never even change in the reference code.

# Progress
## Day 1
- PPM images
- Programmatic image
- Camera + sphere + normal
## Day 2
- Shading: lambertian, metal, glass
- Plane (faked with big sphere)
- World
- Copy world from reference app
- Multisampling
## Day 3
- Progress indicator
- Jimp (for png output)
- Animations
- Multiprocessing

# Running

There are various scripts that can just be run from node directly. For example:

```
node day2-finishbook.js
```

Will output a ppm file that can be viewed in an image viewer or converted with imagemagick:

```
magick convert finishbook.ppm finishbook.png
```

#### Reference Material to keep looking into...:
https://raytracing.github.io/books/RayTracingTheNextWeek.html

https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-ray-tracing/ray-tracing-practical-example.html

https://gabrielgambetta.com/computer-graphics-from-scratch/04-shadows-and-reflections.html

https://www.jakobmaier.at/posts/raytracing/

https://www.kevinbeason.com/smallpt/

https://github.com/noteed/smallpt-hs

<!-- #### vector libs (not used):
https://github.com/maxkueng/victor
https://github.com/tronkko/js-vector
https://github.com/jcoglan/sylvester
https://github.com/toji/gl-matrix -->
