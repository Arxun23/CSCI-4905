import numpy as np
from matplotlib import pyplot as plt
import geopandas as gpd
from pathlib import Path
from multiprocessing import Pool


# loads GPKG file and converts EPSG to 2264
def load_voter_data(data_path: Path) -> gpd.GeoDataFrame:
    print("Loading voter data...")
    pitt_voters = data_path / "pitt_voters.gpkg"
    voter_df = gpd.read_file(pitt_voters)
    return voter_df.to_crs("epsg:2264")


# Filters out voters by parties and UNA voters
def split_voters(voter_df: gpd.GeoDataFrame) -> tuple:
    print("Filtering voter data...")
    dem_voters = voter_df[voter_df["party_cd"] == "DEM"]
    rep_voters = voter_df[voter_df["party_cd"] == "REP"]
    una_voters = voter_df[voter_df["party_cd"] == "UNA"]
    return dem_voters, rep_voters, una_voters


# Saves each voter data into GPKG files
def save_party_voters(data_path: Path, dem_voters: gpd.GeoDataFrame, rep_voters: gpd.GeoDataFrame,
                      una_voters: gpd.GeoDataFrame):
    print("Saving filtered voter data...")
    dem_voters.to_file(data_path / "dem_voters.gpkg", driver="GPKG", OVERWRITE=True)
    rep_voters.to_file(data_path / "rep_voters.gpkg", driver="GPKG", OVERWRITE=True)
    una_voters.to_file(data_path / "UNA_voters.gpkg", driver="GPKG", OVERWRITE=True)


# Process voter data for a certain district
def process_voter_data_for_district(args):
    print("Processing voter data for each district")
    voter_gdf, district_gdf, party, school_level, data_path = args

    # Finds the density
    voter_district_df = counter_voter_by_area(voter_gdf["geometry"], district_gdf["geometry"])
    voter_district_df.insert(0, column="district", value=district_gdf["NAME"])

    # Save to GPKG file
    output_gpkg_file = data_path / f"{party.lower()}_{school_level.lower()}_count_by_district.gpkg"
    voter_district_df.to_file(output_gpkg_file, driver="GPKG")


# Counts voters and finds the density
def counter_voter_by_area(point_series: gpd.GeoSeries, polygons_series: gpd.GeoSeries) -> gpd.GeoDataFrame:

    # makes data frame to save results
    data = np.zeros((len(polygons_series), 2))
    result_df = gpd.GeoDataFrame(data=data, columns=("Count", "Density"), geometry=polygons_series)
    addr_2264 = point_series.to_crs("EPSG:2264")

    # loops through each polygon checking for matches
    for i in polygons_series.index:
        print(f"Processing {i + 1} of {len(polygons_series)} polygons")
        polygon = polygons_series.loc[i]
        for j in point_series.index:
            point = addr_2264.loc[j]
            if point.within(polygon):
                result_df.loc[i, "Count"] += 1.0

        # calculate area
        area_in_sqmi = polygons_series.loc[i].area * 3.86102e-7

        # finds the density
        result_df.loc[i, "Density"] = result_df.loc[i, "Count"] / area_in_sqmi

    return result_df


# finds the density percentage
def analyze_voter_distribution(data_path: Path, school_level: str):

    # Loads voter count data
    dem_file = data_path / f"dem_{school_level.lower()}_count_by_district.gpkg"
    rep_file = data_path / f"rep_{school_level.lower()}_count_by_district.gpkg"
    una_file = data_path / f"una_{school_level.lower()}_count_by_district.gpkg"

    dem_gdf = gpd.read_file(dem_file)
    rep_gdf = gpd.read_file(rep_file)
    una_gdf = gpd.read_file(una_file)

    # finds the percentage
    total_voters_district = dem_gdf["Count"] + rep_gdf["Count"] + una_gdf["Count"]
    dem_per = 100 * dem_gdf["Count"] / total_voters_district
    rep_per = 100 * rep_gdf["Count"] / total_voters_district
    una_per = 100 * una_gdf["Count"] / total_voters_district

    # adds the percentage
    dem_gdf.insert(len(dem_gdf.columns), column="Percent", value=dem_per)
    rep_gdf.insert(len(rep_gdf.columns), column="Percent", value=rep_per)
    una_gdf.insert(len(una_gdf.columns), column="Percent", value=una_per)

    plot_voter_analysis(dem_gdf, rep_gdf, una_gdf, school_level)


# plots out voter results
def plot_voter_analysis(dem_gdf: gpd.GeoDataFrame, rep_gdf: gpd.GeoDataFrame, una_gdf: gpd.GeoDataFrame,
                        school_level: str):

    print(f"Assembling {school_level} school plots")
    # makes a 3x3 plot
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 15))

    # plots voter count for each party
    dem_gdf.plot(column="Count", ax=axes[0, 0], legend=True, cmap="jet")
    rep_gdf.plot(column="Count", ax=axes[1, 0], legend=True, cmap="jet")
    una_gdf.plot(column="Count", ax=axes[2, 0], legend=True, cmap="jet")

    # plots density for each party
    dem_gdf.plot(column="Density", ax=axes[0, 1], legend=True, cmap="jet")
    rep_gdf.plot(column="Density", ax=axes[1, 1], legend=True, cmap="jet")
    una_gdf.plot(column="Density", ax=axes[2, 1], legend=True, cmap="jet")

    # plots percent for each party
    dem_gdf.plot(column="Percent", ax=axes[0, 2], legend=True, cmap="jet")
    rep_gdf.plot(column="Percent", ax=axes[1, 2], legend=True, cmap="jet")
    una_gdf.plot(column="Percent", ax=axes[2, 2], legend=True, cmap="jet")

    # gives a title to each subplot
    axes[0, 0].set_title("DEM Count")
    axes[1, 0].set_title("REP Count")
    axes[2, 0].set_title("UNA Count")

    axes[0, 1].set_title("DEM Density")
    axes[1, 1].set_title("REP Density")
    axes[2, 1].set_title("UNA Density")

    axes[0, 2].set_title("DEM Percentage")
    axes[1, 2].set_title("REP Percentage")
    axes[2, 2].set_title("UNA Percentage")

    # sets the main title
    plt.suptitle("Pitt County Voter Distribution by Precincts, by Luis Sanjuan-Cruz, sanjuancruzl21@students.ecu.edu")

    # saves the figure
    plt.savefig(f"{school_level}_school_distribution_Luis_SanjuanCruz.jpeg")

    # display figure
    print(f"Displaying {school_level} school district charts")
    plt.show()


# Organizes the workflow
def main():
    # Sets up the data paths
    data_path = Path("../Data")
    gis_path = data_path / 'GIS_Data'
    party_group = ("DEM", "REP", "UNA")

    # loads and splits voters and saves them into files
    voter_df = load_voter_data(data_path)
    dem_voters, rep_voters, una_voters = split_voters(voter_df)
    save_party_voters(data_path, dem_voters, rep_voters, una_voters)

    # Create paths for each SHP files for the districts
    districts_files = {
        "elementary": gis_path / "Pitt_County_Elementary_School_Attendance_Districts/"
                                 "Pitt_County_Elementary_School_Attendance_Districts.shp",
        "middle": gis_path / "Pitt_County_Middle_School_Attendance_Districts/"
                             "Pitt_County_Middle_School_Attendance_Districts.shp",
        "high": gis_path / "Pitt_County_High_School_Attendance_Districts/"
                           "Pitt_County_High_School_Attendance_Districts.shp"
    }

    args_list = []

    # prepares arguments for each district
    for school_level, district_file in districts_files.items():
        district_gdf = gpd.read_file(district_file)
        for party in party_group:
            print(f"Processing {party} voters in {school_level} school districts")
            voter_gpkg_file = data_path / f"{party.lower()}_voters.gpkg"
            voter_gdf = gpd.read_file(voter_gpkg_file)
            args_list.append((voter_gdf, district_gdf, party, school_level, data_path))

    # process voter data in parallel
    with Pool() as pool:
        pool.map(process_voter_data_for_district, args_list)

    # calls voter analyzes by school level
    for school_level in ["elementary", "middle", "high"]:
        analyze_voter_distribution(data_path, school_level)


if __name__ == '__main__':
    main()
    print("Done")