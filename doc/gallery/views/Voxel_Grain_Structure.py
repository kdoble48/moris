"""Attached view for the voxel grain structure (3D, nine grains, steady heat conduction).

The default 3D view of a result shows the exterior of the cube colored by temperature, where the
grains (which differ only in conductivity) hardly show. Cut the cube open at mid-depth and color
by element block, so the cut face shows the grain structure the voxel field defines."""


def view(scene):
    scene.clip([0.0, 1.0, 0.0], invert=True)    # keep the half y <= mid-plane; the camera moves to +y
    scene.color_blocks()
    scene.camera(elevation=10, zoom=1.1)
    scene.caption("Cut open at mid-depth and colored by grain (element block) rather than by temperature.")
