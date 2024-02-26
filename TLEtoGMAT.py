from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84
import math
import requests
import re

def mean_to_true_anomaly(M, e, tolerance=1e-10, max_iterations=200):
    """Convert Mean Anomaly to True Anomaly for elliptical orbits."""
    if e == 0:
        # Circular orbit: Mean Anomaly equals True Anomaly
        return M
    # Initial guess for Eccentric Anomaly
    E = M if e < 0.8 else math.pi
    # Newton-Raphson iteration to solve for Eccentric Anomaly
    for _ in range(max_iterations):
        delta = E - e * math.sin(E) - M
        E -= delta / (1 - e * math.cos(E))
        if abs(delta) < tolerance:
            break
    # Convert Eccentric Anomaly to True Anomaly
    true_anomaly = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    return math.degrees(true_anomaly) % 360

def fetch_tle_from_celestrak(url):
    response = requests.get(url)
    response.raise_for_status()  # Raise an error for bad responses
    lines = response.text.strip().split('\n')
    return lines[0], lines[1], lines[2]

def tle_to_gmat(tle_file, gmat_file):
    tle_data = fetch_tle_from_celestrak(satellite_url)
    satellite_name = re.sub(r'\s*\([^)]*\)', '', tle_data[0]).strip()
    line1 = tle_data[1]
    line2 = tle_data[2]
    #with open(tle_file, 'r') as file:
    #    lines = file.readlines()
    #    satellite_name = lines[0].strip()
    #    line1 = lines[1].strip()
    #    line2 = lines[2].strip()

    satellite = twoline2rv(line1, line2, wgs84)
    #print(satellite.no_kozai)    
    # Convert mean motion (revs per day) to radians per minute, then to the SMAprint(satellite.no_kozai)
    mean_motion_rad_min = satellite.no_kozai/60 #* (2 * math.pi / 1440)
    mu = 398600.4418  # km^3/s^2
    sma = (mu / (mean_motion_rad_min ** 2)) ** (1/3)
    
    # Convert Mean Anomaly to True Anomaly
    M = math.radians(satellite.mo)
    e = satellite.ecco
    TA = mean_to_true_anomaly(M, e)
    
    gmat_script = f"""Create Spacecraft {satellite_name};
GMAT {satellite_name}.DateFormat = UTCGregorian;
GMAT {satellite_name}.Epoch = '{satellite.epoch.strftime("%d %b %Y %H:%M:%S.%f")[:-3]}';
GMAT {satellite_name}.CoordinateSystem = EarthMJ2000Eq;
GMAT {satellite_name}.DisplayStateType = Keplerian;
GMAT {satellite_name}.SMA = {sma};
GMAT {satellite_name}.ECC = {e};
GMAT {satellite_name}.INC = {math.degrees(satellite.inclo)};
GMAT {satellite_name}.RAAN = {math.degrees(satellite.nodeo)};
GMAT {satellite_name}.AOP = {math.degrees(satellite.argpo)};
GMAT {satellite_name}.TA = {TA};

Create ForceModel DefaultProp_ForceModel;
GMAT DefaultProp_ForceModel.CentralBody = Earth;
GMAT DefaultProp_ForceModel.PrimaryBodies = {{Earth}};
GMAT DefaultProp_ForceModel.Drag = None;
GMAT DefaultProp_ForceModel.SRP = Off;
GMAT DefaultProp_ForceModel.RelativisticCorrection = Off;
GMAT DefaultProp_ForceModel.ErrorControl = RSSStep;
GMAT DefaultProp_ForceModel.GravityField.Earth.Degree = 4;
GMAT DefaultProp_ForceModel.GravityField.Earth.Order = 4;
GMAT DefaultProp_ForceModel.GravityField.Earth.StmLimit = 100;
GMAT DefaultProp_ForceModel.GravityField.Earth.PotentialFile = 'JGM2.cof';
GMAT DefaultProp_ForceModel.GravityField.Earth.TideModel = 'None';

Create Propagator DefaultProp;
GMAT DefaultProp.FM = DefaultProp_ForceModel;
GMAT DefaultProp.Type = RungeKutta89;
GMAT DefaultProp.InitialStepSize = 60;
GMAT DefaultProp.Accuracy = 9.999999999999999e-12;
GMAT DefaultProp.MinStep = 0.001;
GMAT DefaultProp.MaxStep = 2700;
GMAT DefaultProp.MaxStepAttempts = 50;
GMAT DefaultProp.StopIfAccuracyIsViolated = true;

Create OrbitView DefaultOrbitView;
GMAT DefaultOrbitView.SolverIterations = Current;
GMAT DefaultOrbitView.UpperLeft = [ 0.002941176470588235 0 ];
GMAT DefaultOrbitView.Size = [ 0.5 0.45 ];
GMAT DefaultOrbitView.RelativeZOrder = 20;
GMAT DefaultOrbitView.Maximized = false;
GMAT DefaultOrbitView.Add = {{{satellite_name}, Earth}};
GMAT DefaultOrbitView.CoordinateSystem = EarthMJ2000Eq;
GMAT DefaultOrbitView.DrawObject = [ true true ];
GMAT DefaultOrbitView.DataCollectFrequency = 1;
GMAT DefaultOrbitView.UpdatePlotFrequency = 50;
GMAT DefaultOrbitView.NumPointsToRedraw = 0;
GMAT DefaultOrbitView.ShowPlot = true;
GMAT DefaultOrbitView.MaxPlotPoints = 20000;
GMAT DefaultOrbitView.ShowLabels = true;
GMAT DefaultOrbitView.ViewPointReference = Earth;
GMAT DefaultOrbitView.ViewPointVector = [ 30000 0 0 ];
GMAT DefaultOrbitView.ViewDirection = Earth;
GMAT DefaultOrbitView.ViewScaleFactor = 1;
GMAT DefaultOrbitView.ViewUpCoordinateSystem = EarthMJ2000Eq;
GMAT DefaultOrbitView.ViewUpAxis = Z;
GMAT DefaultOrbitView.EclipticPlane = Off;
GMAT DefaultOrbitView.XYPlane = On;
GMAT DefaultOrbitView.WireFrame = Off;
GMAT DefaultOrbitView.Axes = On;
GMAT DefaultOrbitView.Grid = Off;
GMAT DefaultOrbitView.SunLine = Off;
GMAT DefaultOrbitView.UseInitialView = On;
GMAT DefaultOrbitView.StarCount = 7000;
GMAT DefaultOrbitView.EnableStars = On;
GMAT DefaultOrbitView.EnableConstellations = On;


Create GroundTrackPlot DefaultGroundTrackPlot;
GMAT DefaultGroundTrackPlot.SolverIterations = Current;
GMAT DefaultGroundTrackPlot.UpperLeft = [ 0.002941176470588235 0.4547619047619048 ];
GMAT DefaultGroundTrackPlot.Size = [ 0.5 0.45 ];
GMAT DefaultGroundTrackPlot.RelativeZOrder = 22;
GMAT DefaultGroundTrackPlot.Maximized = false;
GMAT DefaultGroundTrackPlot.Add = {satellite_name};
GMAT DefaultGroundTrackPlot.DataCollectFrequency = 1;
GMAT DefaultGroundTrackPlot.UpdatePlotFrequency = 50;
GMAT DefaultGroundTrackPlot.NumPointsToRedraw = 0;
GMAT DefaultGroundTrackPlot.ShowPlot = true;
GMAT DefaultGroundTrackPlot.MaxPlotPoints = 20000;
GMAT DefaultGroundTrackPlot.CentralBody = Earth;
GMAT DefaultGroundTrackPlot.TextureMap = 'ModifiedBlueMarble.jpg';


BeginMissionSequence;
Propagate DefaultProp({satellite_name}) {{{satellite_name}.ElapsedDays = 1.0}};
"""

    with open(gmat_file, 'w') as file:
        file.write(gmat_script)

# Update the file paths as necessary
# Propagate DefaultProp({satellite_name}) {{}};
# Example URL for the ISS.  Change the "CATNR" number to your satellite's catalog number.
satellite_url = 'https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=TLE'  
gmat_file_path = 'iss_gmat.script'
tle_to_gmat(satellite_url, gmat_file_path)