#
# MIT License
#
# Copyright (c) 2015 Jonathan Shore
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

require (ggplot2)


##
##  Amplitude based labeling of momentum / trend sections
##  ```
##  
##
##  ```
##
AmplitudeBasedLabeler <- R6Class ("AmplitudeBasedLabeler",
    public = list (
        time = NULL,
        prices = NULL,
        cumr = NULL,
        df = NULL,
        minamp = NULL,
        Tinactive = NULL,

        ## constructor
        initialize = function (minamp, Tinactive)
        {
            self$minamp = minamp
            self$Tinactive = Tinactive
        },

        label = function (df, filter = TRUE, time = "time", price = "price")
        {
            self$time <- df[,time]
            self$prices <- df[,price]
            self$cumr <- log(self$prices / self$prices[1]) * 1e4

            # first pass to determine extents
            labels1 <- label_direction (self$cumr, self$minamp, self$Tinactive)
            # refine extents by filtering using a least-squares distance algorithm
            labels2 <- (if (filter)
                filter_direction (self$cumr, labels1, self$minamp)
            else
                labels1)

            self$df <- data.frame (time = self$time, price = self$prices, cumr = self$cumr, label = labels2)
            invisible(self$df)
        },

        ## show plot of prices & label overlay
        plot = function (
            color.up = "#10a4f4",
            color.down = "red",
            color.neutral = "darkgrey",
            size.lines = 0.75,
            dimension = c(10,8),
            title = "")
        {
            data = self$df
            labels = data$label

            ## extract price
            Dprice <- data.frame (time = data$time, value=data$cumr)

            ## segments for 5 diff labels
            up <- data$cumr
            up[labels != +1] <- NA
            Dup <- data.frame (time = data$time, value=up)

            down <- data$cumr
            down[labels != -1] <- NA
            Ddown <- data.frame (time = data$time, value=down)

            neutral <- data$cumr
            neutral[labels != 0] <- NA
            Dneutral <- data.frame (time = data$time, value=neutral)

            ## graph
            options(repr.plot.width=dimension[1], repr.plot.height=dimension[2])
            v <- ggplot() +
                geom_line(aes(x=time, y=value), data=Dprice) +
                geom_line(aes(x=time, y=value), data=Dneutral, color=color.neutral, size=size.lines) +
                geom_line(aes(x=time, y=value), data=Dup, color=color.up, size=size.lines) +
                geom_line(aes(x=time, y=value), data=Ddown, color=color.down, size=size.lines) +
                labs(title=title,y='')

            v
        }
    )
)


