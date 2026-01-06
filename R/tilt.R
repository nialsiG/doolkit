#' Occlusal tilting
#' Account for intra-observer error in orienting the occlusal plane by tilting a mesh along a specified axis.
#' @param mesh object of class mesh3d
#' @param fun a function to apply after each occlusal tilt; the results of this function should be a single numeric value
#' @param axis the axis around which to rotate, either "x", "y" or "z"
#' @param bins the number of non-random rotations (only if randomize = FALSE)
#' @param rotation the maximum tilting range in degrees; tilting will be restricted between rotation and -rotation
#' @param doMap boolean, if TRUE saves a snapshot of a topographic map of the slope for each rotation
#' @param randomize if TRUE randomly tilts within the rotation range
#' @param times the number of random rotations to perform (only if randomize = TRUE)
#' @return From the values calculated after every tilting of the occlusal plane, makes a list with the mean, geometric mean, harmonic mean, variance and standard deviation.
#' @examples
#' mean.slope <- function(mesh) return(mean(slope(mesh)))
#' tilted_slope <- tilt(dkmodel$complex, mean.slope, axis = "x", rotation = 30, times = 100, randomize = FALSE, doMap = FALSE)
#' @export
tilt <- function(mesh, fun, axis = "x", bins = 10, rotation = 30, randomize = FALSE, times = 1000, ...){
  # Perform various checks:
  if (!isa(mesh, what = "mesh3d")) stop("mesh must be an object of class 'mesh3d'")
  if (!is.function(fun)) stop ("fun must be a valid method")

  # Function to get the slope
  WithRotation <- function(mesh, pt1, pt2, theta)
  {
    RotatedMesh <- Morpho::rotaxis3d(mesh, pt1, pt2, theta / 180 * pi)
    FunResult <- fun(RotatedMesh)

    if (length(FunResult) != 1) stop ("fun must return a single value")

    return(round(mean(FunResult), 2))
  }
  # Create list to populate
  Results <- list(average = 0, geommean = 0, harmonicmean = 0, variance = 0, sd = 0, rotations = data.frame(Rotation = 0, Variable = 0))
  # Define the rotation axis
  if (axis == "x")
  {
    pt1 = c(min(mesh$vb[1, ]),
            max(mesh$vb[2, ]) - min(mesh$vb[2, ]),
            max(mesh$vb[3, ]) - min(mesh$vb[3, ]))

    pt2 = c(max(mesh$vb[1, ]),
            max(mesh$vb[2, ]) - min(mesh$vb[2, ]),
            max(mesh$vb[3, ]) - min(mesh$vb[3, ]))
  }
  if (axis == "y")
  {
    pt1 = c(max(mesh$vb[1, ]) - min(mesh$vb[1, ]),
            min(mesh$vb[2, ]),
            max(mesh$vb[3, ]) - min(mesh$vb[3, ]))

    pt2 = c(max(mesh$vb[1, ]) - min(mesh$vb[1, ]),
            max(mesh$vb[2, ]),
            max(mesh$vb[3, ]) - min(mesh$vb[3, ]))
  }
  if (axis == "z")
  {
    pt1 = c(max(mesh$vb[1, ]) - min(mesh$vb[1, ]),
            max(mesh$vb[2, ]) - min(mesh$vb[2, ]),
            min(mesh$vb[3, ]))

    pt2 = c(max(mesh$vb[1, ]) - min(mesh$vb[1, ]),
            max(mesh$vb[2, ]) - min(mesh$vb[2, ]),
            max(mesh$vb[3, ]))
  }
  # Possibility 1: Non-randomized
  if (!randomize)
  {
    #...negative tilting
    for (i in 1:(bins)) {
      theta <- i * (360 / (bins * rotation))
      value <- WithRotation(mesh, pt1, pt2, theta)
      Results$rotations[bins + 1 - i, 1] <- theta
      Results$rotations[bins + 1 - i, 2] <- value
    }
    #...no tilting
    value <- WithRotation(mesh, pt1, pt2, 0)
    Results$rotations[bins + 1, 1] <- 0
    Results$rotations[bins + 1, 2] <- value
    #...positive tilting
    for (j in 1:(bins)) {
      theta <- (- j) * (360 / (bins * rotation))
      value <- WithRotation(mesh, pt1, pt2, theta)
      Results$rotations[(bins + 1 + j), 1] <- theta
      Results$rotations[(bins + 1 + j), 2] <- value
    }
  }
  # Possibility 2: Randomized
  else if (randomize)
  {
    for (i in 1:times)
    {
      theta <- runif(1, min = - rotation, max = rotation + 1)
      value <- WithRotation(mesh, pt1, pt2, theta)
      Results$rotations[i, 1] <- theta
      Results$rotations[i, 2] <- value
    }
  }
  # Global values
  Results$average <- mean(Results$rotations$Variable)
  Results$geommean <- exp(mean(log(Results$rotations$Variable)))
  Results$harmonicmean <- length(Results$rotations$Variable) / sum(1 / Results$rotations$Variable)
  Results$variance <- var(Results$rotations$Variable)
  Results$sd <- sd(Results$rotations$Variable)
  return(Results)
}
