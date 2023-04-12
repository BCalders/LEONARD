# Overview
This code calculates the position of a receiver using Doppler shift measurements from one or more satellites. The receiver's initial position is estimated based on its known coordinates, and the Doppler shift measurements are used to iteratively update the receiver's position estimate.
The develop script is called workingScript
in the Demo folder are a few working implementations with some flaws.

## Initialization
The code begins by defining the initial position of the receiver using its known coordinates. It also initializes a scatter plot to visualize the estimated position of the receiver over time.

## Main Loop
The main loop of the code iteratively updates the receiver's position estimate based on the Doppler shift measurements from the satellite(s). In each iteration of the loop, the following steps are performed:

1. Get the current position and velocity of the satellite(s) using the given satellite states and the current time.

2. Calculate the range vector `rangeVect` between the satellite and the receiver, and normalize it to obtain the unit vector `e`.

3. Calculate the range rate `rhoDot` and range acceleration `rhoDotDot` using the satellite velocity and the unit vector `e`.

4. The expression `(satVel - e .* rhoDot)` represents the relative velocity vector between the satellite and the target. Multiplying this vector by `e` projects it onto the direction of the unit vector pointing from the target to the satellite. This removes any component of the relative velocity vector that is orthogonal to the direction of `e`, leaving only the component of the relative velocity that is parallel to `e`. Dividing this projected velocity by `rho` yields the time derivative of `e`, which is what `edot` represents.

5. Assemble the matrix `H` using the unit vector's derivative `eDot`, and the range acceleration `rhoDotDot`.

6. Calculate the predicted Doppler shift `D_predicted` based on the satellite's frequency and the relative velocity between the satellite and the receiver.

7. Calculate the difference between the observed Doppler shift and the predicted Doppler shift, denoted by `deltaDoppler`.

8. Convert the Doppler shift difference `deltaDoppler` to a range rate difference `deltaD` using the speed of light and the satellite's frequency.

9. Solve the linear equation system `H * deltaX = deltaD` for `deltaX`, the change in the receiver's position estimate.

10. Update the receiver's position estimate by adding `deltaX` to the current estimate.

11. Visualize the updated position estimate on the scatter plot.

The loop continues for a fixed number of iterations or until convergence is achieved.

## Conclusion
Overall, this code implements a simple algorithm for estimating the position of a receiver using Doppler shift measurements from one or more satellites. However, the accuracy and stability of the solution may be limited by the assumptions and simplifications used in the model, as well as the quality of the measurements themselves.