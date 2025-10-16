module DataPrep

using CSV
using DataFrames

function read_asc(filepath)
    # Headers
    headers = [
        "Year",
        "Month",
        "OriginAirportAlpha",
        "OriginAirportCode",
        "OriginWAC",
        "OriginCity",
        "DestinationAirportAlpha",
        "DestinationAirportCode",
        "DestinationWAC",
        "DestinationCity",
        "CarrierAlphaCode",
        "CarrierCode",
        "GroupCode",
        "Distance",
        "ServiceClass",
        "AircraftGroup",
        "AircraftMakeandModel",
        "AircraftConfig",
        "DeparturesPerformed",
        "DeparturesScheduled",
        "AvailableCapacityLbs",
        "AvailableSeats",
        "PassengersTransported",
        "FreightTransported",
        "MailTransported",
        "RampTime",
        "AirbornTime",
        "CarrierWAC",
        "NULL"]

    # Read the ASC file into a DataFrame
    df = CSV.read(filepath, DataFrame, delim='|', header=headers)

    # Get rid of the NULL column
    select!(df, Not(names(df)[end]))

    df = transform(df, :OriginCity => ByRow(s -> split(s, ", ")) => [:OriginCity, :OriginState])
    df = transform(df, :DestinationCity => ByRow(s -> split(s, ", ")) => [:DestinationCity, :DestinationState])


        "CarrierAlphaCode",
        "CarrierCode",
        "GroupCode",
        "Distance",
        "ServiceClass",
        "AircraftGroup",
        "AircraftMakeandModel",
        "AircraftConfig",
        "DeparturesPerformed",
        "DeparturesScheduled",
        "AvailableCapacityLbs",
        "EmptySeats",
        "PassengersTransported",
        "FreightTransported",
        "MailTransported",
        "RampTime",
        "AirbornTime",
        "CarrierWAC",
        "NULL"

    # Define the desired order of columns
    desired_order = [:Year, :Month, :OriginCity, :OriginState, :DestinationCity, :DestinationState,
    :Distance, :RampTime, :AirbornTime,
    :AvailableCapacityLbs, :AvailableSeats, :PassengersTransported, :FreightTransported, :MailTransported,

    :CarrierAlphaCode, :CarrierCode, :GroupCode, :CarrierWAC,

    :AircraftGroup, :AircraftMakeandModel, :AircraftConfig, :ServiceClass,
    :DeparturesPerformed, :DeparturesScheduled,

    :OriginAirportAlpha, :OriginAirportCode, :OriginWAC,
    :DestinationAirportAlpha, :DestinationAirportCode, :DestinationWAC]

    # Create a new DataFrame with reordered columns
    df = df[:, desired_order]

    # Order df by date
    df = sort(df, [:Year, :Month])

    return df
end

function asc_to_csv(filepath)

    data = read_asc(filepath)

    CSV.write("AirplaneDataset.csv", data)

end

function read_coords(filepath, flight_data)
    coords = CSV.read(filepath, DataFrame)

    # get rid of unnecesary columns
    select!(coords, :DISPLAY_AIRPORT_CITY_NAME_FULL, :LATITUDE, :LONGITUDE)
    unique!(coords, :DISPLAY_AIRPORT_CITY_NAME_FULL)

    flight_data.Origin = string.(flight_data.OriginCity,", ", flight_data.OriginState)
    flight_data.Dest = string.(flight_data.DestinationCity,", ", flight_data.DestinationState)

    #TODO: Reduce computation by getting somewhat unique list for flight data

    # Getting coordinates for origins
    origin_coords = innerjoin(flight_data, coords, on = [:Origin => :DISPLAY_AIRPORT_CITY_NAME_FULL], validate = (false, true))
    select!(origin_coords,
        :Origin => :DisplayName,
        :LATITUDE,
        :LONGITUDE,
        :OriginCity => :City,
        :OriginState => :State
    )

    unique!(origin_coords, :DisplayName)

    # Getting coordinates for destinations
    dest_coords = innerjoin(flight_data, coords, on = [:Dest => :DISPLAY_AIRPORT_CITY_NAME_FULL], validate = (false, true))
    select!(dest_coords,
        :Dest => :DisplayName,
        :LATITUDE,
        :LONGITUDE,
        :DestinationCity => :City,
        :DestinationState => :State
    )
    unique!(dest_coords, :DisplayName)

    # Joining the two
    final_coords = vcat(origin_coords, dest_coords)
    unique!(final_coords, :DisplayName)
    return final_coords
end

end # end of DataPrep module