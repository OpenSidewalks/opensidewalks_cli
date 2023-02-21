import json

import numpy as np
from shapely.geometry import mapping


DEFAULT_ZOOM = 19


def write_tasks(tasks, path):
    # Convert into a list of GeoJSON features
    features = []
    for task in tasks:
        # TODO: round the coordinates of the task geometry? We don't need 14
        # decimal places (microscopic) precision and it would decrease the file
        geometry = task["geometry"]

        # Use task 'view' descriptors if available - generate otherwise.
        # x = lon, y = lat, z = map zoom
        x = task.get("x", round(geometry.centroid.x, 6))
        y = task.get("y", round(geometry.centroid.y, 6))
        z = task.get("z", DEFAULT_ZOOM)

        geometry_json = mapping(geometry)
        coords = []
        for ring in geometry_json["coordinates"]:
            new_ring = []
            for coord in ring:
                rounded = list(np.round(coord, 7))
                new_ring.append(rounded)
            coords.append(new_ring)
        geometry_json["coordinates"] = coords

        features.append(
            {
                "type": "Feature",
                "geometry": geometry_json,
                "properties": {"x": x, "y": y, "z": z},
            }
        )

    tasks_fc = {
        "type": "FeatureCollection",
        "features": features,
    }

    with open(path, "w") as f:
        json.dump(tasks_fc, f)
