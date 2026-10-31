import casadi as cs

def batch_fed_bioprocess():
    N = 30

    opti = cs.Opti()
    nx = 7
    nu = 2
    W = opti.variable(nx + nu + 1, N)

    xnp1 = opti.variable(nx, 1)
    tnp1 = opti.variable(1)
    X = cs.hcat([W[0:nx, :], xnp1])
    U = W[nx : nx + nu, :]
    t = cs.hcat([W[nx + nu : nx + nu + 1, :], tnp1])

    # Linear scaling for stability
    #opti.set_linear_scale(t, 10 * cs.DM.ones(t.shape))
    #opti.set_linear_scale(
    #    X,
    #    cs.vcat(
    #        [
    #            cs.DM.ones(1, N + 1),
    #            cs.DM.ones(1, N + 1),
    #            10 * cs.DM.ones(1, N + 1),
    #            cs.DM.ones(1, N + 1),
    #            cs.DM.ones(1, N + 1),
    #            cs.DM.ones(1, N + 1),
    #            cs.DM.ones(1, N + 1),
    #        ]
    #    ),
    #)
    #opti.set_linear_scale(U, cs.DM.ones(U.shape))

    dt = cs.diff(t)
    xs = cs.SX.sym("xs", nx, 1)
    us = cs.SX.sym("us", nu, 1)
    ts = cs.SX.sym("ts")
    t0s = cs.SX.sym("t0s")
    tfs = cs.SX.sym("tfs")

    c1 = 100.0
    c2 = 0.51
    c3 = 4.0

    x1, x2, x3, x4, x5, x6, x7 = (
        xs[0],
        xs[1],
        xs[2],
        xs[3],
        xs[4],
        xs[5],
        xs[6],
    )
    u1, u2 = us[0], us[1]

    t1 = 14.35 + x3 + (x3**2) / 111.5
    t2 = 0.22 + x5
    t3 = (0.22 * x7) / t2 + x6

    g1 = (x3 * t3) / t1
    g2 = (0.233 * x3 * (0.0005 + x5)) / ((0.022 + x5) * t1)
    g3 = (0.09 * x5) / (0.034 + x5)

    dae = {
        "x": xs,
        "t": ts,
        "p": cs.vertcat(us, t0s, tfs),
        "ode": (tfs - t0s)
        * cs.vertcat(
            u1 + u2,
            g1 * x2 - ((u1 + u2) * x2) / x1,
            (u1 / x1) * c1 - ((u1 + u2) * x3) / x1 - (g1 * x2) / c2,
            g2 * x2 - ((u1 + u2) * x4) / x1,
            (u2 * c3) / x1 - ((u1 + u2) * x5) / x1,
            -g3 * x6,
            g3 * (1 - x7),
        ),
    }

    cvodes = cs.integrator(
        "integrator",
        "rk",
        dae,
        0,
        1,
        {"simplify": True, "number_of_finite_elements": 1},
    )

    ires = cvodes.map(N, "serial")(x0=X[:, 0:-1], p=cs.vertcat(U, t[:, 0:-1], t[:, 1:]))
    defect = X[:, 1:] - ires["xf"]

    opti.minimize(-X[3, -1] * X[0, -1])

    opti.subject_to(
        X[:, 0] == cs.DM([1.0, 0.1, 40.0, 0.0, 0.0, 1.0, 0.0])
    )
    opti.subject_to(t[:, 0] == 0)
    opti.subject_to(t[:, -1] == 10.0)

    dtn = (t[:, -1] - t[:, 0]) / N
    opti.subject_to(opti.bounded(0.1 * dtn, dt, 2 * dtn))
    defs=cs.repmat(cs.DM([1.0, 1, 10.0, 1.0, 1.0, 1.0, 1.0]),1,N)
    opti.subject_to(defect == 0)
    opti.subject_to(opti.bounded(0, U, 1))
    opti.subject_to(X[0,:]>=0)
    opti.subject_to(X[4,:]>=0)

    x0=cs.DM([1.0, 0.1, 40.0, 0.0, 0.0, 1.0, 0.0])
    Xg = cs.repmat(x0,1,N+1)
    Ug = cs.vcat([cs.linspace(1,1,N).T]*2)
    tg = cs.linspace(0, 10.0, N + 1).T

    opti.set_initial(t, tg)
    opti.set_initial(X, Xg)
    opti.set_initial(U, Ug)

    opti.solver(
        "ipopt",
    )

    sol = opti.solve()

    from pylab import plot, show, figure, legend, step

    t_val = sol.value(t)
    X_val = sol.value(X)
    U_val = sol.value(U)

    figure(figsize=(10, 6))
    for i in range(nx):
        plot(t_val, X_val[i, :], '-o',label=f"$x_{i+1}$")
    legend()
    figure(figsize=(10, 4))
    step(t_val[:-1], U_val[0, :], label="$u_1$", where="post")
    step(t_val[:-1], U_val[1, :], label="$u_2$", where="post")
    legend()
    show()

    return sol.value(-X[3, -1] * X[0, -1])

batch_fed_bioprocess()