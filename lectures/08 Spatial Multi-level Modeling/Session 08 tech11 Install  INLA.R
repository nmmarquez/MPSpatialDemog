################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 12                #                          
# Install INLA                                                                 #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable", 
                 dep=TRUE)

