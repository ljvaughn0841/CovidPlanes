module DataPrep

using CSV
using DataFrames

function read_route_data(filepath::String)
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

    if endswith(filepath, ".asc")
        # Read the ASC file into a DataFrame
        df = CSV.read(filepath, DataFrame, delim='|', header=headers)
    elseif endswith(filepath, ".csv")
        df = CSV.read(filepath, DataFrame, header=headers)
    else
        println("Unnaceptable Datatype")
        return
    end

    # Get rid of the NULL column
    select!(df, Not(names(df)[end]))

    df = transform(df, :OriginCity => ByRow(s -> split(s, ", ")) => [:OriginCity, :OriginState])
    df = transform(df, :DestinationCity => ByRow(s -> split(s, ", ")) => [:DestinationCity, :DestinationState])

    df.Origin = string.(df.OriginCity,", ", df.OriginState)
    df.Dest = string.(df.DestinationCity,", ", df.DestinationState)

    # Define the desired order of columns
    desired_order = [:Year, :Month, :Origin, :Dest,
    :Distance, :RampTime, :AirbornTime,
    :AvailableCapacityLbs, :AvailableSeats, :PassengersTransported, :FreightTransported, :MailTransported,

    :CarrierAlphaCode, :CarrierCode, :GroupCode, :CarrierWAC,

    :AircraftGroup, :AircraftMakeandModel, :AircraftConfig, :ServiceClass,
    :DeparturesPerformed, :DeparturesScheduled,

    :OriginCity, :OriginState, :DestinationCity, :DestinationState,

    :OriginAirportAlpha, :OriginAirportCode, :OriginWAC,
    :DestinationAirportAlpha, :DestinationAirportCode, :DestinationWAC]

    # Create a new DataFrame with reordered columns
    df = df[:, desired_order]

    # Order df by date
    df = sort(df, [:Year, :Month])

    return df
end

function asc_to_csv(filepath::String)

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



function us_filter(df, US_coords)

    US_coords = filter(:LATITUDE => L -> L > 18.743496240902473, US_coords) # Southern boundary
    US_coords = filter(:LATITUDE => L -> L < 71.55757363763155, US_coords) # Northern Boundary
    US_coords = filter(:LONGITUDE => L -> L > -169.5192533743284, US_coords) # Eastern Boundary
    US_coords = filter(:LONGITUDE => L -> L < -56.99967737291253, US_coords) # Western Boundary

    filter!(:State => S -> S != "PR", US_coords) # Puerto Rico Filter

    US_coords.row_id = 1:nrow(US_coords) # This is necessary for Graphing Points Later

    # Set of valid display names
    valid_names = Set(US_coords.DisplayName)

    # Filter the data DataFrame
    US_df = filter(row -> row.Origin in valid_names && row.Dest in valid_names, df)

    # Removing Self Loops (Flights between areas in the same city)
    filter!(row -> row.Origin != Row.Dest, US_df)

    # Puerto Rico Filter
    filter!(row -> row.Origin != "PR" && row.Dest != "PR", US_df)

    return US_df, US_coords
end





end # end of DataPrep module