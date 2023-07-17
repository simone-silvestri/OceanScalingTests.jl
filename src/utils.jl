using Oceananigans.Units

# Doing the year 1995-2020
leap_year_days(year) = year == 1996 || 
                       year == 2000 || 
                       year == 2004 || 
                       year == 2008 || 
                       year == 2012 || 
                       year == 2016 || 
                       year == 2020 ? UnitRange(1, 29) : UnitRange(1, 28)

monthly_days(year) = [1:31, leap_year_days(year), 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

function realistic_ocean_stop_time(final_year, final_month = 12)
    simulation_days = 0
    # Starts from 01/01/1995
    for year in 1995:final_year-1
        days_in_month = monthly_days(year)
	simulation_days += sum(days_in_month[m][end] for m in 1:12)
    end

    days_in_month = monthly_days(final_year)
    simulation_days += sum(days_in_month[m][end] for m in 1:final_month)

    return simulation_days * days 
end
