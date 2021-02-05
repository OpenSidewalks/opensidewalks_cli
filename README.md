## Fetching / preparing extracts

1. Download OSM data for Ohio from geofabrik: http://download.geofabrik.de/north-america/us/ohio-latest.osm.pbf, place in `data_sources/`

2. Download neighborhoods shapefile from Columbus Open Data Portal: https://opendata.columbus.gov/datasets/columbus-communities, place in `data_sources/`.

3. Open neighborhoods shapefile in QGIS, select 'North Linden' (objectid 383, area name 'North Linden'), save a GeoJSON file with single polygon feature. Save in `intermdiate_data`.

4. Clip openstreetmap data source using neighborhood polygon + the osmium tool: `osmium extract -p intermediate_data/north_linden.geojson data_sources/ohio-latest.osm.pbf -o intermediate_data/north_linden.osm.pbf`

## Creating intersection mapping projects

Run the Python code to extract the intersections and create Voronoi polygons: `poetry run python ./extract_tasks.py /intermediate_data/north_linden.osm.pbf intermediate_data/north_linden.geojson output/north_linden`

## Creating sidewalk mapping projects

