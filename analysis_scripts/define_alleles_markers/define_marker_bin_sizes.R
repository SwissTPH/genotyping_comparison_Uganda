#########################
# Define data frame with bin sizes
# 
# monica.golumbeanu@swisstph.ch
# 24.06.2024
#########################

marker_id = c("K1", "3D7", "FC27", "2490", "MAD20", "R033", "313", 
              "383", "GLURP", "PFPK2", "POLY a", "TA1", "TA109")
bin_size = c(10, 10, 10, 2, 10, 10, 1, 
             1, 50, 2, 2, 2, 2)

marker_bins = cbind.data.frame(marker_id, bin_size)
