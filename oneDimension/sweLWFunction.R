.vFunction <- function(h, u) {
        v1 <- h
        v2 <- h * u
        rbind(v1, v2)
}

.fFunction <- function(h, u, g = 9.81) {
        f1 <- h * u
        f2 <- h * u^2 + 0.5 * g * h^2 
        rbind(f1, f2)
}

.halfStepFunction <- function(h, u, g = 9.81) {
        # function to perform the half-step in LW method
        p <- length(h)
        seqBack <- c(2:p, 1)
        uBackSpace <- u[seqBack]
        hBackSpace <- h[seqBack]
        #
        v <- .vFunction(h = h, u = u)
        vBack <- .vFunction(h = hBackSpace, u = uBackSpace) 
        #
        f <- .fFunction(h = h, u = u, g = g)
        fBack <- .fFunction(h = hBackSpace, u = uBackSpace, g = g)
        #
        vHalfTime <- 0.5 * (vBack + v) - 
                        (dt / 2 * dx) * (fBack - f)
        hOut <- vHalfTime[1, ]
        uOut <- vHalfTime[2, ] / hOut
        list(u = uOut, h = hOut)
}

sweLW <- function(u0 = 0, h0, dt, dx, t.max, g = 9.81) {
        # model 1D swe using FTCS method
        # u0 is initial velocity
        # h0 is initial shape of wave
        # dx is time-step
        # t.max is total time of simulation
        # make sure dt satisfies CFL condition
        if(missing(dt)) dt <- 0.5 * (dx / (max(abs(u0)) + 
                                           sqrt(g * max(abs(h0)))))
        p <- length(h0)
        # total length of wave
        L <- dx * p
        # number of steps
        n <- floor(t.max / dt)
        # initiate
        u <- matrix(ncol = p, nrow = n + 1)
        if (length(u0) == 1) u0 <- rep(u0, times = p)
        u[1, ] <- u0
        h <- matrix(ncol = p, nrow = n + 1)
        h[1, ] <- h0
        seqBack <- c(2:p, 1)
        seqForward <- c(p, 1:(p-1))
        for (i in 1:n) {
                v <- .vFunction(h = h[i, ], u = u[i, ])
                uForwardSpace <- u[i, seqForward]
                hForwardSpace <- h[i, seqForward]
                halfForward <- .halfStepFunction(h = h[i, ],
                                                 u = u[i, ], g = g)
                uHalfForward <- halfForward$u
                hHalfForward <- halfForward$h
                halfBack <- .halfStepFunction(h = hForwardSpace, 
                                              u = uForwardSpace, 
                                              g = g)
                uHalfBack <- halfBack$u
                hHalfBack <- halfBack$h
                fHalfBack <- .fFunction(h = hHalfBack, u = uHalfBack, g = g)
                fHalfForward <- .fFunction(h = hHalfForward, u = uHalfForward,
                                       g = g)
                vForwardTime <- v - (dt / 2 * dx) * (fHalfForward - fHalfBack)
                hTmp <- vForwardTime[1, ]
                uTmp <- vForwardTime[2, ] / hTmp
                u[(i + 1), ] <- uTmp
                h[(i + 1), ] <- hTmp
        }
        list(u = u, h = h)
}

