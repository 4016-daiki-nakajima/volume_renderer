# Volume Renderer

Volume Render with Disney BSDF Implementation on [lajolla renderer](https://github.com/BachiLi/lajolla_public)

![Window](scenes/renders/window.png)
![Window](scenes/renders/disney.png)

# Build
Building is handled by CMake:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

# Run
To run, 
```
./lajolla ../scenes/cbox/cbox.xml
```



## Example Renders
![Teapot cbox](scenes/renders/teapot.png)
![Colored Smoke](scenes/renders/colored_smoke.png)
![Hall](scenes/renders/hall.png)
![Ghosts](scenes/renders/ghosts.png)
![Thanks](scenes/renders/thanks.png)
