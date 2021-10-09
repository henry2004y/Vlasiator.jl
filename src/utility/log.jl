# Runtime log file tracking.

"""
    readlog(file)

Read the run log file, check the iteration status and return the timestamps (exluding the
last) as well as the model running speed in physical seconds per simulated seconds.
"""
function readlog(file)
   isfile(file) || throw(ArgumentError("Cannot open \'$file\': not a file"))

   fid = open(file)
   lines = readlines(fid)
   close(fid)

   nline = findall(startswith("---"), lines)
   steps = lines[nline]
   stepsizes = [parse(Float32, split(step)[10]) for step in steps]
   t = lines[nline .- 1]

   local timestamps
   try
      timestamps = DateTime.(t, dateformat"e u d H:M:S y ")
   catch
      timestamps = DateTime.(t, dateformat"e u  d H:M:S y ")
   end

   timediff = diff(timestamps) / 1000 # duration for each step in [s]
   speed = [timediff[i].value / stepsizes[i] for i in eachindex(timediff)]
   return @views timestamps[1:end-1], speed
end
