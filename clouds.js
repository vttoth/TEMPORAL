// Copyright (c) 2024 Viktor T. Toth
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// A simplistic cloud model

function hash(x, y)
{
  var prime1 = 15731;
  var prime2 = 789221;
  var prime3 = 1376312589;
  var n = x + y * 57;
  n = (n << 13) ^ n;
  return (1.0 - ((n * (n * n * prime1 + prime2) + prime3) & 0x7fffffff) / 1073741824.0);
}

// Linear interpolation function
function lerp(a, b, t)
{
  return a + t * (b - a);
}

// Fade function as defined by Ken Perlin; this eases coordinate values
// so that they will ease towards integral values. This smooths the final output.
function fade(t)
{
  return t * t * t * (t * (t * 6 - 15) + 10);
}

// Gradient function that calculates the dot product between a pseudorandom
// gradient vector and the eight location vectors
function grad(hash, x, y)
{
  var h = hash & 7; // Convert low 3 bits of hash code
  var u = h < 4 ? x : y; // into 8 simple gradient directions,
  var v = h < 4 ? y : x; // and compute the dot product with (x,y).
  return ((h & 1) ? -u : u) + ((h & 2) ? -2 * v : 2 * v);
}

// Perlin noise function
function perlin(x, y)
{
  // Determine grid cell coordinates
  var xi = Math.floor(x) & 255;
  var yi = Math.floor(y) & 255;

  // Relative x, y position in grid cell
  var xf = x - Math.floor(x);
  var yf = y - Math.floor(y);

  // Fade curves for x and y
  var u = fade(xf);
  var v = fade(yf);

  // Hash coordinates of the 4 square corners
  var aa, ab, ba, bb;
  aa = hash(xi, yi);
  ab = hash(xi, yi + 1);
  ba = hash(xi + 1, yi);
  bb = hash(xi + 1, yi + 1);

  // Add blended results from 4 corners of the square
  var x1, x2, y1;
  x1 = lerp(grad(aa, xf, yf), grad(ba, xf - 1, yf), u);
  x2 = lerp(grad(ab, xf, yf - 1), grad(bb, xf - 1, yf - 1), u);
  y1 = lerp(x1, x2, v);

  return (y1 + 1) / 2; // Normalize to [0, 1]
}

// Function to generate fractal noise with multiple octaves
function fractalNoise(x, y, octaves, persistence)
{
  let total = 0;
  let frequency = 1;
  let amplitude = 1;
  let maxValue = 0; // Used for normalizing result to [0, 1]
  for (let i = 0; i < octaves; i++)
  {
    total += perlin(x * frequency, y * frequency) * amplitude;
    maxValue += amplitude;
    amplitude *= persistence;
    frequency *= 2;
  }
  return total / maxValue;
}

// Function to apply contrast to a value
function applyContrast(value, contrastFactor)
{
  // Assumes value is in the range [0, 1]
  return Math.max(0, Math.min(1, contrastFactor * (value - 0.5) + 0.5));
}

// Updated generateCloudCover function with larger clouds
function generateCloudCover(n, threshold)
{
  let cloudCover = new Array(n);
  let baseScale = 1.0 / n; // Base scale for the largest cloud formations
  let detailScale = 2.0 / n; // Detail scale for smaller cloud formations
  let contrastFactor = 2.0; // Adjust this value to control contrast

  for (let i = 0; i < n; i++) {
    cloudCover[i] = new Array(2 * n);
    for (let j = 0; j < 2 * n; j++) {
      // Generate base cloud formations
      let baseNoise = fractalNoise(i * baseScale, j * baseScale, 2, 0.5);
      // Generate smaller cloud details
      let detailNoise = fractalNoise(i * detailScale, j * detailScale, 3, 0.3);
      // Combine the noise layers
      let combinedNoise = (baseNoise + detailNoise) / 2;
      // Apply contrast to the combined noise
      let contrastedNoise = applyContrast(combinedNoise, contrastFactor);
      // Apply a lower threshold to create larger cloud areas
      let cloudValue = contrastedNoise > threshold ? (contrastedNoise - threshold) / (1 - threshold) : 0;
      cloudCover[i][j] = cloudValue;
    }
  }
  return cloudCover;
}

// Function to initialize multiple cloud layers
function initializeCloudLayers(n, layerCount, threshold)
{
  let layers = [];
  for (let i = 0; i < layerCount; i++)
  {
    let direction = Math.random() * Math.PI;
    let layer = {
      cloudCover: generateCloudCover(n, threshold), // Generate the initial cloud cover for this layer
      windSpeed: Math.random() * 0.2, // Random wind speed for this layer
      windDirection: {
//        x: Math.random() * 0.2 - 0.2, // Random wind direction X component
//        y: Math.random() * 0.2 - 0.2  // Random wind direction Y component
         x: Math.cos(direction),
         y: Math.sin(direction)
      }
    };
    layers.push(layer);
  }
  return layers;
}

// Updated progressCloudCover function to move multiple cloud layers
function progressCloudLayers(layers, deltaTime)
{
  const n = layers[0].cloudCover.length;
  layers.forEach(layer => {
    let displacementX = (layer.windDirection.x * layer.windSpeed * deltaTime) % (2 * n);
    let displacementY = (layer.windDirection.y * layer.windSpeed * deltaTime) % n;
    let newCloudCover = new Array(n);

    for (let i = 0; i < n; i++)
    {
      newCloudCover[i] = new Array(2 * n);
      for (let j = 0; j < 2 * n; j++)
      {
        let sourceI = Math.floor((i - displacementY + n) % n);
        let sourceJ = Math.floor((j - displacementX + 2 * n) % (2 * n));
        newCloudCover[i][j] = layer.cloudCover[sourceI][sourceJ];
      }
    }

    layer.cloudCover = newCloudCover; // Update the layer's cloud cover
  });
  return layers;
}

// Function to combine multiple cloud layers into a final cloud cover
function combineCloudLayers(layers, n)
{
  let finalCloudCover = new Array(n).fill(0).map(() => new Array(2 * n).fill(0));

  layers.forEach(layer =>
  {
    for (let i = 0; i < n; i++)
    {
      for (let j = 0; j < 2 * n; j++)
      {
        // Combine the cloud cover values, clamping them to the [0, 1] range
        finalCloudCover[i][j] = Math.min(1, finalCloudCover[i][j] + layer.cloudCover[i][j]);
      }
    }
  });

  return finalCloudCover;
}
