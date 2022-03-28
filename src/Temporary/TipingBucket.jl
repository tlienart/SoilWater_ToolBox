# =============================================================
#		module: tipingBucket
# =============================================================
# module tipingBucket
    using Dates
    import DelimitedFiles
    include("D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\Tool.jl")

    ΔTimeStep = 60.0 * 60.0 # Hourly time step

    Option_BackwardForward = :Forward # <:Backwards> collect ∑Pr between 7h00_8h00; <:Forward> collect cumulative  ∑Prbetween between 8h00_9h00.


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #		FUNCTION : TIPING_BUCKET
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function TIPING_BUCKET()
        Path = "D:\\DATAraw\\LYSIMETERS\\Waikato_data\\ObsPrecipitation\\ObsPrecipitation_Csv\\"

        FileName = ["TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"]

        for iFileName in FileName
            println("===  $iFileName ===")

            iPath_Input = Path * iFileName * "_TipingBucket.csv"

            # READ DATA
				Data = DelimitedFiles.readdlm(iPath_Input, ',')
				Header = Data[1,1:end]
				Data = Data[2:end,begin:end]

				Year, N_Climate = tool.readWrite.READ_HEADER_FAST(Data, Header,"Year")
				Month, ~        = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month")
				Day, ~          = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day")
				Hour, ~         = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour")
				Minute, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute")
				Second, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second")
				Prr, ~           = tool.readWrite.READ_HEADER_FAST(Data, Header, "Rain(mm)")


            # CONVERT INTO DATES
                Date = fill(now()::DateTime,  N_Climate)
                for iT=1:N_Climate
                    Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])
                end # iT=1:N_Climate
            
            # COMPUTE CUMULATIVE PR
                ∑Pr = fill(0.0::Float64, N_Climate)
                ∑Pr[1] = Pr[1]
                for iT=2:N_Climate
                    ∑Pr[iT] = ∑Pr[iT-1] + Pr[iT]
                end # iT=1:N_Climate

            # -------------------------------------------------------------------
            # -------------------------------------------------------------------

            # TIME INTERVAL WHERE WE NEED DATA
                Date_Hourly = []
                Start_Date = DateTime(Year[1], Month[1], Day[1], Hour[1], 0, 0)
        
                push!(Date_Hourly, Start_Date)
                iThourly = 2
                while Date_Hourly[iThourly-1] ≤ Date[N_Climate]
                    push!(Date_Hourly, Date_Hourly[iThourly-1] + Dates.Second(ΔTimeStep))
                    iThourly += 1
                end # iDate ≤ Date[N_Climate]
                N_Hourly = iThourly - 1 

            # GETTING DATA WITH MATCH DATA_HOURLY
                ∑Pr_Hourly = fill(0.0::Float64, N_Hourly)
                iT_Tiping = 1
                for iThourly = 1:N_Hourly
                    while  Date[iT_Tiping] ≤ Date_Hourly[iThourly] && iT_Tiping ≤ N_Climate-1
                        ∑Pr_Hourly[iThourly] = ∑Pr[iT_Tiping]
                        iT_Tiping += 1
                    end # Date_Hourly[iThourly] ≤ Date[iT_Tiping]
                    ∑Pr_Hourly[iThourly] = ∑Pr[iT_Tiping]
                end # iThourly = 1:N_Hourly

            # TRANSFORMING ∑PR TO Pr_Hourly
                Pr_Hourly = fill(0.0, N_Hourly)
                Pr_Hourly[1] = ∑Pr_Hourly[1]

                for iThourly = 2:N_Hourly
                    Pr_Hourly[iThourly] = ∑Pr_Hourly[iThourly] - ∑Pr_Hourly[iThourly-1]
                end

            # IF OPTION FORWARD : SHIFTING DATA 
                if Option_BackwardForward == :Forward
                    Date_Hourly = Date_Hourly[1:end-1]
                    Pr_Hourly = Pr_Hourly[2:end]
                end 

            # SAVING TO FILE
                iPath_Output = Path * "Output\\" * iFileName * "_TipingBucket_Output.csv"

                Header = ["Year", "Month", "Day", "Hour", "Minute", "Second", "Pr[mm]"]

                open(iPath_Output, "w") do io
                    DelimitedFiles.writedlm(io,[Header] , ",",) # Header
                    DelimitedFiles.writedlm(io, [year.(Date_Hourly) month.(Date_Hourly) day.(Date_Hourly) hour.(Date_Hourly) minute.(Date_Hourly) second.(Date_Hourly) Pr_Hourly], ",")
                end # open
        end # for iFileName in FileName

    return nothing
    end  # function: TIPING_BUCKET
    
# end  # module: tipingBucket
# ............................................................

TIPING_BUCKET()