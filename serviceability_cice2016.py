import numpy as np
from beam import getrcbeam_cice2016, getfrpbeam_cice2016, getfrpcreep_cice2016

# load bridge girder
rcbeam = getrcbeam_cice2016()
frpbeam = getfrpbeam_cice2016()
frpcreep = getfrpcreep_cice2016()
beamarray = np.array([rcbeam,frpbeam,frpcreep])

for ib,beamEle in enumerate(beamarray):
    # basic data
    bf = 2600.; hf = 190.; b = 400.; h = 790.; L=9100.
    fc = beamEle.mat.fc
    frein = beamEle.mat.fr
    Mdl = 222.5e6
    # Mll = 550.4e6
    Mll = 289.63e6    #AASHTO: three design trucks divided equally among 5 girders
    [Rm, c] = beamEle.getmoment()
    if ib==0:
        Arein = np.sum(beamEle.geo.As)
        d = np.average(beamEle.geo.xs,weights=beamEle.geo.As)
        Erein = 200e3
        Ma = Mll
        # Ma = Mdl+Mll
    elif ib==1:
        Arein = np.sum(beamEle.geo.Afb)
        d = np.average(beamEle.geo.xf,weights=beamEle.geo.Afb)
        Erein = 45e3
        # Ma = Mdl+Mll
        Ma = Mll
    elif ib==2:
        Arein = np.sum(beamEle.geo.Afb)
        d = np.average(beamEle.geo.xf,weights=beamEle.geo.Afb)
        Erein = 45e3
        Ma = Mdl
    ecu = 0.003
    eru = frein/Erein
    Ec = 4700.*np.sqrt(fc)

    fr = 0.62*1.*np.sqrt(fc)
    cy = h - (h**2*b+hf**2*(bf-b))/(2*(bf*hf+(h-hf)*b))
    Ig = 1./12.*b*h**3+b*h*(cy-h/2)**2+1./12*(bf-b)*hf**3+(bf-b)*hf*(h-cy-hf/2)**2
    yt = h/2    #conservatively
    Mcr = fr*Ig/yt

    rhof = Arein/(b*d)
    nf = Erein/Ec
    Icr = 1./3*bf*c**3+nf*np.sum(Arein)*(d-c)**2

    if ib == 1:
        beta1 = 0.85
        rhofb = 0.85*beta1*fc/frein*(Erein*ecu)/(Erein*ecu+frein)
        betad = 1./5*(rhof/rhofb)
        if betad>1: betad =1.
    else:
        betad = 1.
    if Mcr<=Ma:
        Ie = (Mcr/Ma)**3*betad*Ig + (1-(Mcr/Ma)**3)*Icr
    else:
        Ie = (Mcr/Ma)**3*betad*Ig
    if Ie>Ig: Ie=Ig
    df = 5*Ma*L**2/(48*Ec*Ie)

    if ib == 0:
        print 'RC span-to-deflection = {}'.format(L/df)
    else:
        print 'FRP span-to-deflection = {}'.format(L/df)
