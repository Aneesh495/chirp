# ADR-0004: Doppler-Aware Ellipsoidal Metric for Point Cloud Clustering

## Status
Accepted

## Context
Standard 3D LiDAR point cloud clustering relies entirely on Cartesian spatial proximity (Euclidean distance $d = \sqrt{\Delta x^2 + \Delta y^2 + \Delta z^2}$). In automotive radar, however, spatial resolution is significantly coarser than LiDAR, causing distinct obstacles in adjacent traffic lanes or oncoming vehicles to appear spatially adjacent or overlapping in range and azimuth.

A critical failure mode of Euclidean DBSCAN on radar point clouds is under-segmentation: two vehicles traveling in opposite directions at 50 km/h and -50 km/h that pass within 2 meters of each other are merged into a single monolithic bounding box. This creates dangerous kinematic track state instability, throwing Kalman filter velocity estimates into chaotic oscillations.

Radar measurements provide an independent, high-precision physical dimension that LiDAR lacks: instantaneous radial Doppler velocity ($v_r$). Detections belonging to the same rigid body must share consistent Doppler velocities, while detections from adjacent obstacles traveling at different speeds exhibit massive velocity separation.

## Decision
We implement a 4D Doppler-Aware Density-Based Spatial Clustering of Applications with Noise (Velocity-Aware DBSCAN) engine that groups points based on a coupled spatio-Doppler ellipsoidal distance metric:

$$D(p_1, p_2) = \sqrt{\frac{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}{\epsilon_{\text{range}}^2} + \frac{(v_1 - v_2)^2}{\epsilon_{\text{velocity}}^2}} \le 1.0$$

1. Coordinates are normalized into a dimensionless 4D metric space: $\tilde{p} = [x/\epsilon_r, y/\epsilon_r, z/\epsilon_r, v/\epsilon_v]^T$.
2. Spatial indexing and neighbor discovery are accelerated using a $k$-d tree (`SpatialKdTree`) operating on normalized 4D points.
3. Queue expansion utilizes a strictly ordered FIFO traversal with hash/flag membership verification to guarantee deterministic, reproducible clustering across all execution environments.

## Consequences
### Positive
- Completely eliminates under-segmentation of adjacent vehicles traveling at distinct velocities.
- Correctly separates crossing pedestrians, cyclists, and moving vehicles even when their physical positions overlap in 3D Cartesian coordinates.
- Greatly stabilizes downstream multi-target track filtering (IMM-UKF), yielding tight innovation bounds and reducing track loss rates.

### Negative
- Rotational motion of long vehicles (e.g. articulated trucks turning an intersection) produces Doppler gradients across the vehicle body, which requires setting $\epsilon_{\text{velocity}}$ sufficiently wide ($\approx 2.0\text{ to }2.5\text{ m/s}$) to avoid over-segmenting single physical vehicles.
