STREET_TAGS = [
    "primary",
    "secondary",
    "tertiary",
    "residential",
    "unclassified",
    "trunk",
]


# Loop through the OSM file and extract network + coordinate information
def street_filter(d):
    # Exclude non-street/pedestrian ways
    if "highway" not in d:
        return False

    # Include main street types
    if d["highway"] in STREET_TAGS:
        return True

    # Include alleys
    if d["highway"] == "service":
        if "service" not in d:
            # Is a 'bare' service road, might have a crossing.
            return True
        elif d["service"] == "alley":
            # Is an alley
            return True
        elif d["service"] == "drive-through":
            # Include drive-thrus
            return True

    # Exclude everything else
    return False
