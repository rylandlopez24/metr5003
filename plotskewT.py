def plotskewT(temp,pres):
    import matplotlib.pyplot as plt
    g = 9.81 # m/s^2
    Cp = 1004 # J/kgK 
    Rd = 287 # J/kgK
    Lv = 2.5e6 # J/kg
    Rv = 461.5 # J/kgK

    fig, ax = plt.subplots(figsize=(xsize, ysize))


    ax.semilogy() #log based axis for pressure
    ax.set_ylim(1000, 100)
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter()) #getting both out of scientific notation
    ax.set_ylabel('Pressure (hPa)')

    ax.set_xlim(-40, 50)
    ax.set_xlabel('Temperature (°C)')




    for i in range(-110, 50, 10): # every 10 degc
        ax.plot(skewT(np.array([i, i]), np.array([1000, 100])), np.array([1000, 100]), color='black', linewidth=0.5, linestyle='--')




    for thetaq in range(-30, 240, 20):
        ax.plot(skewT(isentrope(thetaq), np.arange(1000, 100, -1)), np.arange(1000, 100, -1), color='red', linewidth=0.25, label='Dry Adiabats' if thetaq == -30 else None) # only label the first one to avoid clutter


    return fig, ax

def isentrope(theta0):

    g = 9.81 # m/s^2
    Cp = 1004 # J/kgK 
    Rd = 287 # J/kgK
    Lv = 2.5e6 # J/kg
    Rv = 461.5 # J/kgK

    theta = theta0 + 273.15 # converting to kelvin
    p0 = 1000
    pf = 100
    oneisentrope = np.arange(p0, pf, -1) # basically making for one single isentrope, the pressure values from 1000 to 100 hPa and finding 
    # temp given pressure and potential temp 

    return theta * (oneisentrope/p0)**(Rd/Cp) - 273.15 # basically solving for temperature from theta

def moist(T): # for an initial temperature
    g = 9.81 # m/s^2
    Cp = 1004 # J/kgK 
    Rd = 287 # J/kgK
    Lv = 2.5e6 # J/kg
    Rv = 461.5 # J/kgK
    p0 = 1000 # hPa
    Tk = T + 273.15 # converting to kelvin
    pressure = p0

    finaltemps = []

    while pressure > 100:

        es = 6.112 * np.exp(17.67 * (Tk - 273.15)/(Tk - 29.65)) # saturation vapor pressure in hPa

        rs = (Rd/Rv) * (es/(pressure - es)) # saturation mixing ratio in kg/kg

        gamma = (Rd*Tk/(Cp*pressure))*((1+(Lv*rs/(Rd*Tk)))/(1+(Lv**2*rs/(Cp*Rv*Tk**2)))) # moist adiabatic lapse rate in K/hPa - adjusted it by using hydrostatic equation to get rid of dz
        finaltemps.append(Tk - 273.15) # appending the temperature in celsius to the list
        Tk -= gamma # updating the temperature in kelvin
        pressure -= 1
    return finaltemps

def skewT(temp, pres, skew=-40): #I like 40 - defining skew here
    return (temp) + skew * np.log(pres/1000.0)
