// Model parameters

var WIDTH = (new URLSearchParams(window.location.search)).get('width');
if (!WIDTH) WIDTH = 360;                 // Source map pixel width
var HEIGHT = WIDTH / 2;
var RADIUS = HEIGHT / 2;
var INTERVAL = 50;                       // ms for animation
var EPOCH = 2460086;                     // May 21, 2023 at 12:00:00
var ME = {LAT: 45, LON: -75, J: EPOCH};  // Direction to observer
var IMSIZE = 1300;                       // Image size in units of aperture
{
  let meLat = (new URLSearchParams(window.location.search)).get('lat');
  let meLon = (new URLSearchParams(window.location.search)).get('lon');
  let meJ = (new URLSearchParams(window.location.search)).get('J');
  if (meLat) ME.LAT = meLat;
  if (meLon) ME.LON = meLon;
  if (meJ) ME.J = meJ;
}

// Constants

let SIDAY = 23.93446972222222;  // Sidereal day in hours

// Global objects

let data = [];         // The map data, cylindrical projection read upon window load
let image = [];        // The projected globe
let blur = [];         // The blurred globe
let sol = [];          // The deblurred solution

var observations = []; // "Observations" of the blurred globe will be collected here

// Canvas objects
var cia = null;        // The "CIA map" cylindrical projection showing light and shadows
var globe = null;      // The projected globe before blurring
var bglobe = null;     // The blurred globe
var sglobe = null;     // The deblurred, recovered globe

var timer = 0;         // Global timer for running animations


// Initializing blank arrays
for (let i = 0; i < HEIGHT; i++) data[i] = Array(WIDTH).fill(0);
for (let i = 0; i < HEIGHT; i++) image[i] = Array(HEIGHT).fill(0);
for (let i = 0; i < HEIGHT; i++) blur[i] = Array(HEIGHT).fill(0);

// Helpers and physics functions

// The averaged PSF (elliptic integral) approximation
function E(r)
{
  //r *= 1300/HEIGHT;     // Actual pixel size for a ~30pc planet
  r *= IMSIZE/HEIGHT;     // Actual pixel size for a ~1pc planet
  let r2 = r/2;
  return ((1-Math.exp(-Math.PI*r2))/Math.PI+Math.exp(-(r*r)/Math.PI))/(1+r2);
}

// A sigmoid-type helper to make images appear more visually appealing
function sigmoid(z)
{
  z = 4*z-2;
  return 1 / (1 + Math.exp(-z));
}

// Return the azimuth and elevation of the Sun at LAT, LON, epoch J
function getIllumination(LAT, LON, J)
{
  const J2000 = 2451545.0;  // Julian date of the J2000 epoch
  const MSL = 280.46;       // Mean Solar Longitude, degrees
  const MDM = 0.9856474;    // Mean Daily Motion of the Sun along the ecliptic, in days
  const MSA = 357.528;      // Mean Solar Anomaly, degrees
  const MDE = 0.9856003;    // Mean Daily motion of the Sun in the Ecliptic, in days
  const AMP = 1.915;        // Amplitude of largest term equation of the center, measuring Earth's orbital ellipticity
  const AM2 = 0.020;        // Amplitude of second largest term of the above
  const MOB = 23.439;       // Mean obliquity of the ecliptic, degrees
  const ROB = 0.0000004;    // Rate of change of obliquity, degrees per day
  const GMT = 100.46;       // Greenich Mean Sidereal Time at J2000, in degrees

  const n = J-J2000;
  const L = MSL+MDM*n;
  const g = (MSA+MDE*n)*Math.PI/180;
  const l = (L+AMP*Math.sin(g)+AM2*Math.sin(2*g))*Math.PI/180;
  const e = (MOB-ROB*n)*Math.PI/180;
  const a = Math.atan2(Math.cos(e)*Math.sin(l),Math.cos(l));
  const d = Math.asin(Math.sin(e)*Math.sin(l));
  const UT = (J-Math.floor(J)-0.5)*24;
  const LST = (GMT+MDM*n+LON+15*UT)*Math.PI/180;
  const H = LST-a;
  const A = 180/Math.PI*Math.atan2(Math.sin(H),Math.cos(H)*Math.sin(LAT*Math.PI/180)-Math.tan(d)*Math.cos(LAT*Math.PI/180));
  const h = 180/Math.PI*Math.asin(Math.sin(LAT*Math.PI/180)*Math.sin(d)+Math.cos(LAT*Math.PI/180)*Math.cos(d)*Math.cos(H));

  return {A: A, h: h}; // Solar azimuth, solar elevation
}

// The main map

// Read the map data from a file and display
function readMap()
{
  const request = new XMLHttpRequest();
  request.open('GET', 'marble' + WIDTH + 'x' + HEIGHT + '.raw');
  request.responseType = 'arraybuffer';

  request.onload = () =>
  {
    const arrayBuffer = request.response;
    const buffer = new Float64Array(arrayBuffer);
  
    let index = 0;
    for (let i = 0; i < HEIGHT; i++)
    {
      for (let j = 0; j < WIDTH; j++)
      {
        data[i][j] = buffer[index++];
      }
    }

    const canvas = document.createElement('canvas');
    canvas.width = WIDTH;
    canvas.height = HEIGHT;

    const ctx = canvas.getContext('2d');
    let img = ctx.createImageData(WIDTH, HEIGHT);
    let pixels = img.data;

    for (let i = 0; i < HEIGHT; i++)
    {
      for (let j = 0; j < WIDTH; j++)
      {
        const value = data[i][j];
        let p = Math.round(value*255);
        pixels[i*WIDTH*4 + j*4 + 0] = p;
        pixels[i*WIDTH*4 + j*4 + 1] = p;
        pixels[i*WIDTH*4 + j*4 + 2] = p;
        pixels[i*WIDTH*4 + j*4 + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);

    // Append canvas to document
    document.getElementById('orig').appendChild(canvas);
    canvas.style.width = 720;
    canvas.style.height = 360;
    canvas.style.marginRight = '20px';

    showIt(ME.LAT, ME.LON, (ME.J - EPOCH)*24);
  };

  request.send();
}

// The projected globe

// Calculate the globe
function doMap(f0, l0, h)
{
  f0 = -f0;             // Because y goes from top to bottom in the data set ...

  J = EPOCH + h / 24.0; // For illumination

  h *= 24/SIDAY;        // Convert to sidereal 'hours' (1/24th a full rotation)
  h = (h%24+24)%24-24;  // To ensure that we are in the proper range

  for (let i = -RADIUS; i < RADIUS; i++)
  {
    for (let j = -RADIUS; j < RADIUS; j++)
    {
      let r = Math.sqrt(i**2 + j**2);
      if (r > RADIUS) { image[i+RADIUS][j+RADIUS] = 0; continue; }

      let d = Math.asin(r/RADIUS);
      let t = -Math.atan2(j, i); //Math.PI/2 - Math.atan2(i, j);

       // =      ASIN(     SIN(A2*  PI()/180)*COS(C2/6371)+COS(A2*   PI()/180)*    SIN(C2/6371)*COS(D2*PI()/180))*180/PI()
      let f = Math.asin(Math.sin(f0*Math.PI/180)*Math.cos(d)+Math.cos(f0*Math.PI/180)*Math.sin(d)*Math.cos(t));
       // = B2+     180/PI()   *  ATAN2(COS(C2/6371)-SIN(A2*   PI()/180)*SIN(E2*PI()/180),SIN(D2*PI()/180)*SIN(C2/6371)*COS(A2*PI()/180))
      let l = l0-15*h+180/Math.PI*Math.atan2(Math.cos(d)-Math.sin(f0*Math.PI/180)*Math.sin(f),Math.sin(t)*Math.sin(d)*Math.cos(f0*Math.PI/180));

      //let p = (Math.floor(f*HEIGHT/Math.PI) + RADIUS) % HEIGHT;
      //let q = (Math.floor(RADIUS+l*WIDTH/360)+WIDTH) % WIDTH;
      //let H = (getIllumination(90-p*180/HEIGHT,q*360/WIDTH-180,J).h*Math.PI/180);
      let p = (f*180/Math.PI + 90) % 180;
      let q = (l + 450) % 360;
      let H = (getIllumination(90-p,q-180,J).h*Math.PI/180);
      p = Math.floor(p*HEIGHT/180);
      q = Math.floor(q*WIDTH/360);
      let I = H>0 ? Math.sin(H) : 0;
      image[i+RADIUS][j+RADIUS] = I * (0.25 + 0.75*data[p][q]);
    }
  }
}

// Show the globe
function showMap()
{
  if (!globe)
  {
    globe = document.createElement('canvas');
    globe.width = HEIGHT;
    globe.height = HEIGHT;

    // Append canvas to document
    document.getElementById('spinner').appendChild(globe);
    globe.style.width = 720;
    globe.style.height = 720;
  }

  const ctx = globe.getContext('2d');

  let img = ctx.createImageData(HEIGHT, HEIGHT);
  let pixels = img.data;

  for (let i = 0; i < HEIGHT; i++)
  {
    for (let j = 0; j < HEIGHT; j++)
    {
      let value = image[i][j];
      value = 1 - (1 - value)**2;
      let p = Math.round(value*255);
      pixels[i*HEIGHT*4 + j*4 + 0] = p;
      pixels[i*HEIGHT*4 + j*4 + 1] = p;
      pixels[i*HEIGHT*4 + j*4 + 2] = p;
      pixels[i*HEIGHT*4 + j*4 + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
}

// Show the "CIA map"
function showCIA(J)
{
  if (!cia)
  {
    cia = document.createElement('canvas');
    cia.width = WIDTH;
    cia.height = HEIGHT;
    cia.style.width = 720;
    cia.style.height = 360;
    cia.style.marginRight = '20px';

    document.getElementById('cia').appendChild(cia);
  }

  const ctx = cia.getContext('2d');

  let img = ctx.createImageData(WIDTH, HEIGHT);
  let pixels = img.data;

  for (let i = 0; i < HEIGHT; i++)
  {
    for (let j = 0; j < WIDTH; j++)
    {
    const value = sigmoid(getIllumination(90-i*180/HEIGHT,j*360/WIDTH-180,J).h*Math.PI/180) * (0.25 + 0.75*data[i][j]);

    let p = Math.round(value*255);
    pixels[i*WIDTH*4 + j*4 + 0] = p;
    pixels[i*WIDTH*4 + j*4 + 1] = p;
    pixels[i*WIDTH*4 + j*4 + 2] = p;
    pixels[i*WIDTH*4 + j*4 + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
}

// The blurred globe

// Run the blur in a worker thread
function doBlur(e)
{
  context = e.data;
  var IMSIZE = context.IMSIZE;
  var HEIGHT = context.HEIGHT;
  var blur = context.blur;
  var image = context.image;
  var E = new Function ('IMSIZE', 'HEIGHT', 'r', context.E + "; return E(r);");

  if (!self.isRunning)
  {
    for (let i = 0; i < HEIGHT; i++)
    {
      for (let j = 0; j < HEIGHT; j++)
      {
        blur[i][j] = 0;
        for (k = 0; k < HEIGHT; k++)
        {
          for (l = 0; l < HEIGHT; l++)
          {
            let r = Math.sqrt((i-k)**2 + (j-l)**2);
            blur[i][j] += E(IMSIZE, HEIGHT, r) * image[k][l];
          }
        }
      }
      self.postMessage({progress: i/HEIGHT});
    }
    self.postMessage(blur);
  }
}

// Initialize and run the worker thread, manage UI
function runBlur()
{
  if (document.isBlurring) return;
  document.isBlurring = true;
  document.body.style.cursor = 'wait';
  document.querySelectorAll('.controls button').forEach((e) => {e.disabled = true; e.style.cursor = 'inherit'; });

  var worker = new Worker(window.URL.createObjectURL(new Blob(["onmessage=" + doBlur.toString()], {type: "text/javascript"})));
  worker.onmessage = function(e)
  {
    var data = e.data;
    if (data.progress != null)
    {
      progress = (100*data.progress).toFixed(1);
      document.getElementById('progress').innerHTML = "<progress value='" + progress + "' max='100'></progress>&nbsp;" + progress + "%";
    }
    else
    {
      blur = data;
      document.getElementById('progress').innerHTML = '';
      showBlur();
      document.isBlurring = false;
      document.body.style.cursor = '';
      document.querySelectorAll('.controls button').forEach((e) => {e.disabled = false; e.style.cursor = ''; });
      bglobe.disabled = false;
      bglobe.style.cursor = 'crosshair';
      worker.terminate();
    }
  };

  worker.postMessage({IMSIZE: IMSIZE, HEIGHT: HEIGHT, blur: blur, image: image, E: E.toString()});
}

// Display the blurred globe
function showBlur()
{
    if (!bglobe)
    {
      bglobe = document.createElement('canvas');
      bglobe.width = HEIGHT;
      bglobe.height = HEIGHT;

      document.getElementById('blurred').appendChild(bglobe);
      bglobe.style.width = 720;
      bglobe.style.height = 720;

      bglobe.addEventListener("click", onBlurClick);
    }

    const bctx = bglobe.getContext('2d');

    let img = bctx.createImageData(HEIGHT, HEIGHT);
    let pixels = img.data;
    for (let i = 0; i < HEIGHT; i++)
    {
      for (let j = 0; j < HEIGHT; j++)
      {
        //let value = blur[i][j] / HEIGHT;
        //value = 1 - (1 - value)**2;
        //let p = Math.round(value*255);
		let p = blur[i][j]*IMSIZE/HEIGHT;
        pixels[i*HEIGHT*4 + j*4 + 0] = p;
        pixels[i*HEIGHT*4 + j*4 + 1] = p;
        pixels[i*HEIGHT*4 + j*4 + 2] = p;
        pixels[i*HEIGHT*4 + j*4 + 3] = 255;
      }
    }
    bctx.putImageData(img, 0, 0);
}

// Manage click events, collect "observations"
function onBlurClick(e)
{
  if (bglobe.disabled) return;

  let i = Math.floor(e.offsetY * HEIGHT / parseInt(bglobe.style.height));
  let j = Math.floor(e.offsetX * HEIGHT / parseInt(bglobe.style.height));
  console.log("bglobe[" + i + "][" + j + "]");
  let found = false;  // To avoid duplicates/degenerate convolution matrix
  observations.every((o) =>
  {
    if (o.i == i && o.j == j && Math.abs(o.J - ME.J) < 1e-4)
    {
      found = true;
      return false;
    }
    return true;
  });
  if (!found)
  {
    observations.push({J: ME.J, i: i, j: j, blur: blur[i][j]});
    document.getElementById("nobs").innerText = observations.length;
  }
}

// Events, animations

// Skip forward or backward -- blur must be redone
function skip(delta)
{
  ME.J += delta/24;
  showIt(ME.LAT, ME.LON, (ME.J - EPOCH)*24);
  if (bglobe)
  {
    bglobe.disabled = true;
    bglobe.style.cursor = 'not-allowed';
  }
}

// Compute and show the globe given direction and time
function showIt(lat, lon, h)
{
  doMap(lat, lon, h);
  showMap();
  showCIA(EPOCH + h/24);
}

// Animate the globe (not the blur) -- does not update ME
function runIt(lat, lon, epoch, dur, step = 1)
{
  if (timer) return;
  var h = 0;
  timer = setInterval(function()
  {
    showIt(lat, lon, epoch + h);
    h += step;
    if (h > dur)
    {
      clearInterval(timer);
      timer = 0;
    }
  }, INTERVAL);
}

// Animate one sidereal day at a time
function runSD(lat, lon, epoch, dur)
{
  runIt(lat, lon, epoch, dur, SIDAY);
}

// Initialize map on window load
window.addEventListener('DOMContentLoaded', () =>
{
  readMap();
  document.getElementById('counter').value = (stepCount.toString()).replace(/\s+/g,' ');
  document.getElementById('divisor').value = (stepDivisor.toString()).replace(/\s+/g,' ');
});

function setCounter()
{
  eval("stepCount = " + document.getElementById('counter').value);
}

function setDivisor()
{
  eval("stepDivisor = " + document.getElementById('divisor').value);
}

// Load and save

function saveAll()
{
  observations.forEach((o) => { delete o.S; });
  let savedata = JSON.stringify(
  {
      WIDTH: WIDTH,
      ME: ME,
      observations: observations,
      counter: stepCount.toString(),
      divisor: stepDivisor.toString()
  });
  let blob = new Blob([savedata], {type: "application/json"});
  let url = URL.createObjectURL(blob);
  let a = document.createElement("a");
  a.href = url;
  a.download=document.title + ".json";
  a.click();
}

function loadAll()
{
  let fileInput = document.createElement("input");
  fileInput.type = "file";
  fileInput.accept = ".json, application/json";
  fileInput.addEventListener("change", () =>
  {
    let file = fileInput.files[0];
    fileName = file.name.replace(/.json$/,"");
    let reader = new FileReader();
    reader.addEventListener("load", () =>
    {
      let data = JSON.parse(reader.result);
      WIDTH = data.WIDTH;
      ME = data.ME;
      observations = data.observations;

      if (data.counter)
      {
        document.getElementById('counter').value = data.counter;
        setCounter();
      }
      if (data.divisor)
      {
        document.getElementById('divisor').value = data.divisor;
        setDivisor();
      }

      document.getElementById("nobs").innerText = observations.length;
      fileInput.remove();
    });
    reader.readAsText(file);
  });
  fileInput.click();
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Deblur -- experimental code

function showDeblur(X, H, W)
{
  if (!sglobe)
  {
    sglobe = document.createElement('canvas');
    sglobe.style.width = 720;
    sglobe.style.height = 270;
    document.getElementById('deblur').appendChild(sglobe);
  }
  sglobe.width = W;
  sglobe.height = H;

  var sctx = sglobe.getContext('2d');

  let img = sctx.createImageData(W, H);
  let pixels = img.data;
  let base = Math.min(...X);
  let scale = Math.max(...X) - base;
//  let base = X.reduce((a, b) => a + b, 0) / X.length;
//  let scale = Math.sqrt(X.reduce((a, b) => a + (b - base)**2, 0) / X.length);
//  base -= 2*scale;
//  scale *= 4;

  for (let i = 0; i < H; i++)
  {
    for (let j = 0; j < W; j++)
    {
      let value = X[i*W + j];
      value = (value - base) / scale;
//if (value<0) value=0;
//if (value>1) value=1;
      let p = Math.round(value*255);
      pixels[i*W*4 + j*4 + 0] = p;
      pixels[i*W*4 + j*4 + 1] = p;
      pixels[i*W*4 + j*4 + 2] = p;
      pixels[i*W*4 + j*4 + 3] = 255;
    }
  }
  sctx.putImageData(img, 0, 0);
}

// Deblur fitting worker thread
function workFit(e)
{
  var context = e.data;
  var ME = context.ME;
  var SIDAY = context.SIDAY;
  var STEPDIV = context.STEPDIV;
  var KMAX = context.KMAX;
  var RADIUS = context.RADIUS;
  var EPOCH = context.EPOCH;
  var observations = context.observations;
  var IMSIZE = context.IMSIZE;
  var HEIGHT = context.HEIGHT;
  var E = new Function ('IMSIZE', 'HEIGHT', 'r', context.E + "; return E(r);");
  var getIllumination = new Function('LAT', 'LON', 'J', context.getIllumination + "; return getIllumination(LAT, LON, J);");

  // Compute the right size for the deconvolved pixel array
  function getSize(N)
  {
    let solution;
    while (!solution && N > 0)
    {
      // Find all factor pairs that satisfy aspect ratio
      let factors = [];
      for (let h = 1; h <= N; h++)
      {
        if (N % h === 0)
        {
          let w = N / h;
          if (h / w >= 1/4 && h / w <= 1/2)
          {
            factors.push([h, w]);
          }
        }
      }
  
      // Calculate squareness metric for each factor pair
      let squareness = factors.map(f => Math.min(...f) / Math.max(...f));
  
      // Find index of factor pair with maximum squareness
      let maxIndex = squareness.indexOf(Math.max(...squareness));
  
      // If tie, choose factor pair using most data points
      if (squareness.filter(s => s === squareness[maxIndex]).length > 1)
      {
        let maxArea = 0;
        for (let i = 0; i < factors.length; i++)
        {
          if (squareness[i] === squareness[maxIndex])
          {
            if (factors[i][0] * factors[i][1] > maxArea)
            {
              maxIndex = i;
              maxArea = factors[i][0] * factors[i][1];
            }
          }
        }
      }
      solution = factors[maxIndex];
      N--;
    }
    return solution;
  }

  // Run the model that determines the mapping
  function model()
  {
    var f0 = -ME.LAT;
    var l0 = ME.LON;
    let p0 = (ME.LAT >= 0) ? 0 : -ME.LAT;

    // We create a (slightly) overdetermined system on purpose
    let theSize = getSize(Math.floor(observations.length * 0.95));
    if (!theSize) return;

    let H = theSize[0];
    let W = theSize[1];

    console.log("Model height: " + H + ", width: " + W);

    let count = 0;
    let total = observations.length;
    observations.forEach((o) =>
    {
      self.postMessage({progress: (count++)/(2*total)});
      o.S = Array(H*W).fill(0);

      let h = (o.J - EPOCH) * 24;
      h *= 24/SIDAY;        // Convert to sidereal 'hours' (1/24th a full rotation)
      h = (h%24+24)%24-24;  // To ensure that we are in the proper range

      for (let i = -RADIUS; i < RADIUS; i++)
      {
        for (let j = -RADIUS; j < RADIUS; j++)
        {
          let r = Math.sqrt((o.i-i-RADIUS)**2 + (o.j-j-RADIUS)**2);
          let blur = E(IMSIZE, HEIGHT, r);

          r = Math.sqrt(i**2 + j**2);
          if (r > RADIUS) continue;
          let d = Math.asin(r/RADIUS);
          let t = -Math.atan2(j, i); //Math.PI/2 - Math.atan2(i, j);
          let f = Math.asin(Math.sin(f0*Math.PI/180)*Math.cos(d)+Math.cos(f0*Math.PI/180)*Math.sin(d)*Math.cos(t));
          let l = l0-15*h+180/Math.PI*Math.atan2(Math.cos(d)-Math.sin(f0*Math.PI/180)*Math.sin(f),Math.sin(t)*Math.sin(d)*Math.cos(f0*Math.PI/180));

          let p = (f * 180 / Math.PI + 90) % 180;
          let q = (l + 450) % 360;

          // Promote by half a pixel width and height to find the center
          p += (180 - Math.abs(ME.LAT)) / (2*H);
          q += 360 / (2*W);
          let I = Math.sin(getIllumination(90-p,q-180,o.J).h*Math.PI/180);
          if (I <= 0) continue;
          // Make sure we subtract the half pixel by which positions were promoted
          q = Math.floor(q * W / 360 - 0.5);
          p = Math.floor((p - p0) * H / (180 - Math.abs(ME.LAT)) - 0.5);
          let k = p * W + q;
          o.S[k] += I * blur;
        }
      }
    });

    return {W: W, H: H};
  }

  // Solve for the surface
  function constrainedLeastSquares(C, V, A, B)
  {
    let M = C.length; 
    let N = C[0].length;
    let X = [];
  
    // Initialize X to midpoint of range
    for (let i = 0; i < N; i++) { X[i] = (A + B) / 2; }
  
    // Calculate residuals
    let R = [];
    for (let i = 0; i < M; i++)
    {
      let sum = 0;
      for (let j = 0; j < N; j++) { sum += C[i][j] * X[j]; }
      R[i] = V[i] - sum;
    }

    // Update X using steepest descent
    let delta = 1e99;
    let olddelta;
    let K = 0;
    do
    {
      olddelta = delta;
      delta = 0;
      for (let j = 0; j < N; j++)
      {
        let grad = 0;
        for (let i = 0; i < M; i++) { grad += C[i][j] * R[i]; }
        let step = grad / STEPDIV;
        X[j] = Math.min(Math.max(X[j] + step, A), B);
        delta += Math.abs(step);
      }
    
      // Recalculate residuals
      for (let i = 0; i < M; i++)
      {
        let sum = 0;
        for (let j = 0; j < N; j++) { sum += C[i][j] * X[j]; }
        R[i] = V[i] - sum;
      }
      self.postMessage({delta: delta, progress: (KMAX+K)/(2*KMAX)});
    } while (delta > 1e-6 && delta < olddelta && ++K < KMAX);
    console.log("" + (K) + " iterations, delta=" + delta);

    return X;
  }

  self.postMessage({progress: 0});

  var m = model();
  sol = constrainedLeastSquares(observations.map(o => o.S),observations.map(o => o.blur),-600,600);

  self.postMessage({H:m.H, W:m.W, sol: sol});

}

// These helpers are to be revised
function stepDivisor()
{
  return observations.length;
}

function stepCount()
{
  return 100;
}

// Used for evaluating the goodness of the fit
function resample(data, newWidth, newHeight)
{
  let width = data.length;
  let height = data[0].length;
  
  let xRatio = width / newWidth; 
  let yRatio = height / newHeight;
  
  let newData = [];
  for (let y = 0; y < newHeight; y++)
  {
    let newRow = [];
    for (let x = 0; x < newWidth; x++)
    {
      let x1 = Math.floor(x * xRatio);
      let x2 = Math.ceil(x * xRatio);
      let y1 = Math.floor(y * yRatio);
      let y2 = Math.ceil(y * yRatio);
      
      let topLeft = data[x1][y1];
      let topRight = data[x2][y1];
      let bottomLeft = data[x1][y2];
      let bottomRight = data[x2][y2];
      
      let xLerp = (xRatio * x) - x1;
      let yLerp = (yRatio * y) - y1;
      let top = topLeft + (topRight - topLeft) * xLerp;
      let bottom = bottomLeft + (bottomRight - bottomLeft) * xLerp;
      let newVal = top + (bottom - top) * yLerp;
      
      newRow.push(newVal);
    }
    newData.push(newRow);
  }
  return newData;
}

function doFit()
{
//  var m = model();
//  sol = constrainedLeastSquares(observations.map(o => o.S),observations.map(o => o.blur),-600,600);
//  showDeblur(sol, m.H, m.W);

  if (document.isFitting) return;
  document.isFitting = true;
  document.body.style.cursor = 'wait';
  document.querySelectorAll('.controls button').forEach((e) => {e.disabled = true; e.style.cursor = 'inherit'; });

  var delta;

  var worker = new Worker(window.URL.createObjectURL(new Blob(["onmessage=" + workFit.toString()], {type: "text/javascript"})));
  worker.onmessage = function(e)
  {
    var result = e.data;
    if (result.progress != null)
    {
      delta = result.delta;
      progress = (100*result.progress).toFixed(1);
      document.getElementById('progress').innerHTML = "<progress value='" + progress + "' max='100'></progress>&nbsp;" + progress + "%";
    }
    else
    {
      sol = result.sol;
      document.getElementById('progress').innerHTML = '';
      showDeblur(sol, result.H, result.W);
      document.isFitting = false;
      document.body.style.cursor = '';
      document.querySelectorAll('.controls button').forEach((e) => {e.disabled = false; e.style.cursor = ''; });
      worker.terminate();

      let error = (resample(data,sglobe.width,sglobe.height).flat().map((d,i) => (d - sol[i])**2)).reduce((a,b)=>a+b)/sol.length;

      document.getElementById('comment').innerHTML = "Model height: " + sglobe.height + ", width: " + sglobe.width + "<br/>" +
                                                     "Residual: " + delta.toPrecision(3) + ", Error: " + error.toPrecision(3);
    }
  }

  worker.postMessage({ME: ME, SIDAY:SIDAY, STEPDIV: stepDivisor(), KMAX: stepCount(), RADIUS:RADIUS, EPOCH: EPOCH, observations: observations,
                      IMSIZE: IMSIZE, HEIGHT: HEIGHT, E: E.toString(), getIllumination : getIllumination.toString()});

}
