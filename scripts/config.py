# config.py  —  YOUR settings. Claude does NOT overwrite this file, so your
# token and monitor list are safe when the code (poller.py) gets updated.

# 1) bearer token (already filled in). Expires ~Aug 26; when it dies, log into
#    dash.aqi.in, grab a fresh token from any request's authorization header,
#    and replace the string below.
TOKEN = "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJodHRwOi8vYWlycXVhbGl0eS5hcWkuaW4vYXBpL3YxL1NpZ25JbldpdGhHb29nbGUiLCJpYXQiOjE3ODc1NzIyNDgsImV4cCI6MTc4ODc4MTg0OCwibmJmIjoxNzg3NTcyMjQ4LCJqdGkiOiJRdFRnRnVOUWF1N24wQWpsIiwic3ViIjoiMTM3OTc1IiwicHJ2IjoiMjNiZDVjODk0OWY2MDBhZGIzOWU3MDFjNDAwODcyZGI3YTU5NzZmNyJ9.66nDEyYvcByFEtv13zaKZMUIB3jzjysQTvERylynWwo"

# 2) PM2.5 parameter id (SAME for every device). sensorid 2 = AQI-US, so find the
#    PM2.5 one: in the dashboard click a device's PM2.5 metric with DevTools open,
#    find the GetDeviceLastHourData request whose response says "sensorName":"PM2.5",
#    and read its sensorid header. Put that number here.
PM25_SENSORID = "3"

# 2b) other parameters to also pull (each is a sensorid, same for all devices).
#     Find each the same way: click that metric in the dashboard, read the
#     sensorid header of the request whose response names it. Leave blank to skip.
PARAMS = {
    "pm25":        "3",
    "pm1":         "5",
    "pm10":        "4",
    "aqi":         "2",
    "temperature": "11",
    "humidity":    "12",
}

# 3) how often to poll, in seconds (60 = every minute, 600 = every 10 min)
INTERVAL_SEC = 60

# 4) one serialno per monitor. Click each monitor in the dashboard, read the
#    serialno header, and add it here mapped to your label. Set the anchor's role
#    to "anchor" and any not-yet-deployed ones to "pending" (those are skipped).
MONITORS = {
    "1":  {"serialno": "3844BE864090", "role": "anchor"},   # anchor, by MIDC
    "2":  {"serialno": "3844BE867328", "role": "field"},     # online
    "3":  {"serialno": "3844BE854D6C", "role": "field"},
    "4":  {"serialno": "3844BE854E54", "role": "field"},     # polled; poller auto-skips it while offline
    "5":  {"serialno": "3844BE8582E4", "role": "field"},
    "6":  {"serialno": "3844BE854044", "role": "field"},
    "7":  {"serialno": "3844BE8540B8", "role": "field"},     # polled; poller auto-skips it while offline
    "8":  {"serialno": "441BF6035300", "role": "field"},     # env board dead (temp 146C, RH stuck 100): PM kept, temp/RH synthesized from neighbors
    "9":  {"serialno": "3844BE83181C", "role": "field"},
    "10": {"serialno": "3844BE85898C", "role": "field"},
}

# 5) per-unit calibration.
#    corrected_pm25 = b0 + b_raw*raw  (+ b_rh*RH  once a reference-fit RH term exists)
#    Monitor-to-monitor harmonization fit from the Aug 3-10 Civil Lines colocation:
#    each unit regressed onto the network-median PM2.5 (despiked). r >= 0.996 for all
#    units, RMSE <= ~1 ug/m3. This removes inter-unit bias (relative calibration).
#    It does NOT fix the shared high-RH monsoon inflation - that needs an absolute
#    reference (gravimetric filters), added later as a "b_rh" key per unit.
CALIB = {
    "1":  {"b0":  0.114, "b_raw": 1.0026},
    "2":  {"b0":  0.231, "b_raw": 0.9832},
    "3":  {"b0": -0.246, "b_raw": 1.0233},
    "4":  {"b0":  0.132, "b_raw": 0.9851},
    "5":  {"b0":  0.313, "b_raw": 0.9242},
    "6":  {"b0":  0.243, "b_raw": 0.9593},
    "7":  {"b0":  0.289, "b_raw": 0.9907},
    "8":  {"b0":  0.308, "b_raw": 1.0085},
    "9":  {"b0":  0.165, "b_raw": 1.0216},
    "10": {"b0":  0.489, "b_raw": 0.9847},
}

# 6) INTERIM humidity correction - kappa-Kohler hygroscopic growth.
#    Optical sensors read high in humid air (water condenses on particles). This
#    divides that growth back out:
#        PM_dry = PM / (1 + RH_KAPPA * aw/(1-aw)),   aw = min(RH, RH_CAP)/100
#    It is a PHYSICAL, literature-based correction, NOT validated against a local
#    reference (the MIDC reference is dead; gravimetric filters not yet weighed), so
#    treat it as provisional. Tune the knob below:
#      - RH_KAPPA = 0 disables it entirely (back to monitor-to-monitor only).
#      - Lower it toward ~0.2 if the map looks over-corrected; Prana/Plantower units
#        already compensate for humidity internally, so full physical kappa can be too
#        strong. Raise toward 0.5 for a more hygroscopic aerosol.
#    Replace with a reference-fit b_rh term (in CALIB) once the filters come back.
RH_CORRECTION = True
RH_KAPPA = 0.0      # OFF for now. Rooftop siting isolates the plant signal (not local
                    # dust/traffic), so a large downward humidity correction isn't
                    # warranted, and it's unvalidated until the gravimetric filters exist.
                    # Over-reporting is the safe direction for a public health map, so we
                    # keep monitor-to-monitor only. Once the filters are weighed, either
                    # set a reference-fit "b_rh" per unit in CALIB, or raise this kappa.
RH_CAP   = 92       # %; cap RH in the growth factor to avoid near-saturation blow-up

# 7) OpenAQ reference node - a live coverage node pulled from the OpenAQ v3 API.
#    Chauhan Colony, Chandrapur (CPCB/MPCB reference-grade). It is NOT run through the
#    monitor-to-monitor or RH calibration - it is already a reference instrument, so its
#    value passes through as-is. Data is ~15-min with a 1-2 h lag. Set OPENAQ_KEY="" to
#    disable. Add more reference stations to OPENAQ_REFS as {id, location_id, lat, lon}.
OPENAQ_KEY  = "fef4863aeec1745fbfcace0b20a19c5c4e25269bed550704d824a0a3d69fb9e9"
OPENAQ_REFS = [
    {"id": "R1", "location_id": 3409317, "lat": 19.9629, "lon": 79.298714,
     "name": "Chauhan Colony (CPCB ref)"},
]
