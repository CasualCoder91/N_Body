import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy import units as u

PISKUNOV = ascii.read("Piskunov_2007-OC_masses.txt")


def plot_radius_distribution():

    core_rad = PISKUNOV["rc"]
    tidal_rad = PISKUNOV["rt"]

    plt.figure()
    plot_series_histogram(core_rad, 30, label="core", alpha=0.5)
    plot_series_histogram(tidal_rad, 30, label="tidal", alpha=0.5)
    plt.legend(loc=1)

    plt.xlabel("Cluster radii [pc]")
    plt.title("Open cluster radii distribution (Piskunov+07)")
    plt.savefig("young_cluster_radii.png")


def plot_mass_distribution():
    mass = 10 ** PISKUNOV["logM"]

    plt.figure()
    plot_series_histogram(mass, 30)

    plt.xlabel("Cluster mass [M$_\odot$]")
    plt.title("Open cluster mass distribution (Piskunov+07)")
    plt.savefig("young_cluster_masses.png")


def plot_age_distribution():
    age = 10 ** PISKUNOV["logt"]

    plt.figure()
    plot_series_histogram(age/1e6, 30)

    plt.xlabel("Cluster age [Myr]")
    plt.title("Open cluster age distribution (Piskunov+07)")
    plt.savefig("young_cluster_ages.png")


def plot_distance_distribution():
    distance = PISKUNOV["Dist"]

    plt.figure()
    plot_series_histogram(distance, 30)

    plt.xlabel("Cluster distance [pc]")
    plt.title("Open cluster distance distribution (Piskunov+07)")
    plt.savefig("young_cluster_distances.png")


def plot_coordinates():
    glat = PISKUNOV["GLAT"]
    glon = PISKUNOV["GLON"]
    dist = PISKUNOV["Dist"]

    from astropy.coordinates import SkyCoord
    coords = SkyCoord(glon*u.deg, glat*u.deg, dist*u.pc, frame="galactic")
    plt.plot(coords.cartesian.x/1e3, coords.cartesian.y/1e3, ".")
    plt.scatter([0], [0], s=100, c="r")
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.xlabel("X [kpc]")
    plt.ylabel("Y [kpc]")
    plt.gca().set_aspect("equal")
    plt.savefig("young_cluster_coordinates.png")


def plot_series_histogram(x, bins, logs="x", **kwargs):

    if "x" in logs:
        logx = np.log10(x)
        log_mu = np.average(logx)
        mu = 10**log_mu
        log_sig = np.std(logx)
        sigs = 10**np.array([log_mu-log_sig, log_mu+log_sig])
        if isinstance(bins, int):
            bins = np.logspace(np.min(logx), np.max(logx), bins)
    else:
        mu = np.average(x)
        sig = np.std(x)
        sigs = 10**np.array([mu-sig, mu+sig])
        if isinstance(bins, int):
            bins = np.linspace(np.min(x), np.max(x), bins)

    plt.hist(x, bins=bins, **kwargs)

    alpha = 0.7
    for i, s in enumerate(sigs):
        text = f"{i*2-1}$\sigma$ = {np.round(s, 2)}"
        plt.axvline(s, c="k", ls=":", zorder=100, alpha=alpha)
        plt.text(s, 0, text, rotation=90, horizontalalignment="right",
                 verticalalignment="bottom", color="k", alpha=alpha)

    text = f"$\mu$ = {np.round(mu, 2)}"
    plt.axvline(mu, c="k", ls="--", zorder=101, alpha=alpha)
    plt.text(mu, 0, text, rotation=90, horizontalalignment="right",
             verticalalignment="bottom", color="k", alpha=alpha)

    if "x" in logs:
        plt.semilogx()


# plot_mass_distribution()
# plot_radius_distribution()
# plot_age_distribution()
# plot_distance_distribution()
plot_coordinates()