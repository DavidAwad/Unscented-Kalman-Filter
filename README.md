# Unscented Kalman Filter

This project is a simple working imlementation of an unscented kalman filter, it is able to make accurate estimations about an objects position and velocity and adjust it's predictions based on additional data that it processes and updates.


The UKF is founded on the intuition that it is easier to approximate a probability distribution that it is to approximate an arbitrary nonlinear function or transformation.
The sigma points are chosen so that their mean and covariance to be exactly x a k−1 and Pk−1.
Each sigma point is then propagated through the nonlinearity yielding in the end a cloud of transformed points.
The new estimated mean and covariance are then computed based on their statistics. This process is called unscented transformation.
The unscented transformation is a method for calculating the statistics of a random variable which undergoes a nonlinear transformation


---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/f437b8b0-f2d8-43b0-9662-72ac4e4029c1)
for instructions and the project rubric.
