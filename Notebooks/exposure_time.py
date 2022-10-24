import numpy as np
from scipy import integrate

# hard coded variables
h = 6.626*10**-34 # planck constant
c = 2.998*10**8 # speed of light in a vacuum
M_sun = 4.74 # absolute magnitude of the Sun
L_sun = 3.828*10**26 # luminosity of the Sun
au = 1.496*10**11 # astronomical unit
pixel_size = 3.8*10**-6 # pixel size of camera
max_photons = 32000 # maximum number of photons allowed to land on a pixel
f_sun = L_sun/(4*np.pi*au**2) # flux of the Sun
m_sun = M_sun-5+5*np.log10(au/(3.086*10**16)) # apparent magnitude of the Sun
fwhm_arcseconds = 4 # full width at half maximum of starlight incident on the camera, spread as a Gaussian

# function for gaussian distribution
def gauss(x, std):
    return (1/(std*np.sqrt(2*np.pi)))*np.exp((-1/2)*((x)/std)**2)

def exposure_time(m_star, telescope, filter):

    # if statements for the different telescopes. the telescopes have different diameters and focal lengths
    if telescope == 'Apo refractor':
        diameter = 130 * 10 ** -3
        focal_length = 910 * 10 ** -3
    elif telescope == 'Ritchey-Chretien':
        diameter = 200 * 10 ** -3
        focal_length = 1600 * 10 ** -3
    else:
        print("Choose one of the two available telescopes, Apo refractor or Ritchey-Chretien.")

    # if statements for the type of filter, which determines the wavelength of photons incident on the camera, and the percentage of photons transmitted
    if filter == 'L':
        T = 0.95
        l1 = 390 * 10 ** -9
        l2 = 700 * 10 ** -9
    elif filter == 'R':
        T = 0.96
        l1 = 590 * 10 ** -9
        l2 = 700 * 10 ** -9
    elif filter == 'G':
        T = 0.95
        l1 = 490 * 10 ** -9
        l2 = 570 * 10 ** -9
    elif filter == 'B':
        T = 0.93
        l1 = 390 * 10 ** -9
        l2 = 500 * 10 ** -9
    elif filter == 'Halpha':
        T = 1
        l1 = 656.3 * 10 ** -9
        l2 = 656.3 * 10 ** -9
    elif filter == 'OIII':
        T = 1
        l1 = 500 * 10 ** -9
        l2 = 500 * 10 ** -9
    elif filter == 'SII':
        T = 1
        l1 = 672 * 10 ** -9
        l2 = 672 * 10 ** -9
    elif filter == 'None':
        T = 1
        l1 = 551 * 10 ** -9
        l2 = 551 * 10 ** -9
    else:
        print("Choose a filter, L, R, G, B, Halpha, OIII, SII or NO for no filter")

    # determining the width of the Gaussian spread of light
    angular_resolution = (pixel_size / focal_length) * (180 / np.pi) * 3600  # angular resolution of the camera (arcseconds per pixel)
    fwhm_pixels = fwhm_arcseconds / angular_resolution  # number of pixels the full width at half maximum covers

    # using the full width at half maximum to calculate the standard deviation of the starlight entering the camera
    std = fwhm_pixels / (np.sqrt(8 * np.log(2)))

    # area of telescope aperture
    area = np.pi * (diameter / 2) ** 2

    # calculation of the flux of the star
    f_star = f_sun * 100 ** ((m_sun - m_star) / 5)

    # approximating the average energy of a photon let through by the choice of filter
    wavelength = np.linspace(l1, l2, 10000)  # range of wavelengths spanning the full width at half maximum
    E_photon = h * c / wavelength  # energy of a photon
    total_E = sum(E_photon)  # energy of different photon wavelengths summed up
    ave_E = total_E / 10000  # average energy of a photon

    # determining the number of photons incident on the camera
    f_photon = f_star / ave_E  # photon flux
    n_photons = area * f_photon * T  # number of photons that enter the camera

    # calculating the greatest probability of a photon reaching a specific pixel using a gaussian distribution for incident photons
    p1d = integrate.quad(gauss, -1 / 2, 1 / 2, args = (std,))  # integrating the gaussian function
    p2d = p1d[0] ** 2  # greatest probability of a photon hitting a specific pixel

    # final calculations for the exposure time for a specific apparent magnitude
    pixel_photons = n_photons * p2d  # number of photons landing on the pixel per second
    t = max_photons / pixel_photons  # exposure time
    tm = t / 60  # exposure time in minutes
    th = tm / 60  # exposure time in hours
    hours = int(th)
    minutes = int(tm - 60 * hours)
    seconds = round(t - 60 * minutes - 3600 * hours)
    print(t)
    print('the exposure time needed for 32000 photons to hit one pixel is', round(t), 'second(s) or', hours, 'hour(s),',
          minutes, 'minute(s)', seconds, 'and second(s)')

exposure_time(10.7, 'Ritchey-Chretien', 'L') # input apparent magnitude of star, telescope and filter here
