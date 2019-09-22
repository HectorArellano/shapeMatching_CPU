import {Particle} from './Particle.js';

//Generate the corresponding particles
let particles = [];
let container = document.querySelector('.container');
let radius = 8;

for(let j = 0; j < 128; j ++) {
    for (let i = 0; i < 128; i++) {
        //if( Math.abs(i - 64) < 10 && Math.abs(j - 64) < 10) {
        let r = Math.pow((i - 64), 2) + Math.pow((j - 64), 2);
        if(r < 32 * 32 && i < 50) {
            particles.push(new Particle(container, i, j, particles.length));
            particles[particles.length - 1].radius(radius);
        }
    }
}

let restDensity = 1000;

//Mass is calculated based on m = ρ * v where v(m3) = 1m * 1m * (1 / 128)m and
//The total mass is divided by the amount of particles allocated in the space.

let particleMass = restDensity;
let searchRadius = 1.5;
let relaxParameter = 0.095;
let deltaTime = 0.03;

let wConstant = (315 / (64 * Math.PI * Math.pow(searchRadius, 9)));
let densityConstant = wConstant * particleMass;
let gradWconstant = -45 / (Math.PI * Math.pow(searchRadius, 6));

let tensileConstant = 40;
let tensilePower = 4;
let tensileDistance = 0.1 * searchRadius;

let center = {x:64, y:64};
let circleRadius = 64 * 0.494;

let searchSquare = searchRadius * searchRadius;
let iterations = 10;
let looping = true;
let viscosityConstant = 0.1 * 45 / (Math.PI * Math.pow(searchRadius, 6) * restDensity);

let voxels = [];
let totalVoxels = 128 * 128;
for (let i = 0; i < totalVoxels; i++) voxels[i] = [];
let maxParticlesPerVoxel = 4;


console.log("total particles: " + particles.length);
console.log("particle mass: " + particleMass);
console.log("poly6 constant: " + wConstant);
console.log("spiky gradient constant: " + gradWconstant);


let iterativeStep = () => {

    //Evaluate constrains for each particle
    for(let particle_i of particles) {

        let sum_k_grad_Ci = 0.;
        let grad_pi_Ci = {x:0, y:0};
        particle_i.density = 0;
        particle_i.neighbors = [];

        for (let u = -1; u <= 1; u++) {
            for (let w = -1; w <= 1; w++) {

                let index = Math.floor(particle_i.x) + u + 128 * (Math.floor(particle_i.y) + w);

                for(let v = 0; v < Math.min(voxels[index].length, maxParticlesPerVoxel); v ++) {

                    let particle_j = particles[voxels[index][v]];

                    let distance = particle_i.distance(particle_j);

                    if(distance < searchRadius) {

                        particle_i.density += Math.pow(searchSquare - distance * distance, 3);

                        if(distance > 0) {

                            let partial = searchRadius - distance;
                            let vector = {x: (particle_i.x - particle_j.x), y: (particle_i.y - particle_j.y)}
                            let multiplier = gradWconstant * partial * partial / (restDensity * distance);

                            vector.x *= multiplier;
                            vector.y *= multiplier;

                            grad_pi_Ci.x += vector.x;
                            grad_pi_Ci.y += vector.y;
                            sum_k_grad_Ci += vector.x * vector.x + vector.y + vector.y;
                        }
                    }
                }
            }
        }



        sum_k_grad_Ci += grad_pi_Ci.x * grad_pi_Ci.x + grad_pi_Ci.y * grad_pi_Ci.y;

        particle_i.density *= densityConstant;
        particle_i.constrain = particle_i.density / restDensity - 1;
        particle_i.lambda = - particle_i.constrain / (sum_k_grad_Ci + relaxParameter);
    }

    //Evaluate displacements for each particle
    for(let particle_i of particles) {

        let deltaDisplacement_i = {x:0, y:0}

        for (let u = -1; u <= 1; u++) {
            for (let w = -1; w <= 1; w++) {

                let index = Math.floor(particle_i.x) + u + 128 * (Math.floor(particle_i.y) + w);

                for (let v = 0; v < Math.min(voxels[index].length, maxParticlesPerVoxel); v++) {

                    let particle_j = particles[voxels[index][v]];

                    let distance = particle_i.distance(particle_j);

                    if (distance > 0 && distance < searchRadius) {

                        //Lambda correction
                        let lambdaCorrection = (searchSquare - distance * distance) / (searchSquare - tensileDistance * tensileDistance);
                        lambdaCorrection = -tensileConstant * Math.pow(lambdaCorrection, 3 * tensilePower);

                        let partial = searchRadius - distance;
                        let vector = {x: (particle_i.x - particle_j.x), y: (particle_i.y - particle_j.y)};

                        let m = ((particle_i.lambda + particle_j.lambda + lambdaCorrection) * partial * partial) / distance;

                        deltaDisplacement_i.x += vector.x * m;
                        deltaDisplacement_i.y += vector.y * m;
                    }
                }
            }
        }

        particle_i.x += gradWconstant * deltaDisplacement_i.x / restDensity;
        particle_i.y += gradWconstant * deltaDisplacement_i.y / restDensity;


        let normal = {x: particle_i.x - center.x, y: particle_i.y - center.y}
        let normalL = Math.sqrt( normal.x * normal.x + normal.y * normal.y);
        var dist = normalL - circleRadius;
        var lMax = 0.0001;
        if(dist > 0.) {
            normal.x /= normalL;
            normal.y /= normalL;
            var cP = {x: center.x + circleRadius * normal.x, y: center.y + circleRadius * normal.y}
            particle_i.x = cP.x - lMax * normal.x;
            particle_i.y = cP.y - lMax * normal.y;
        }

        particle_i.displacement = deltaDisplacement_i;

    }
}

let simulationStep = () => {

    if(looping) requestAnimationFrame(simulationStep);

    //Set the particles inside the voxels
    for (let i = 0; i < totalVoxels; i++) voxels[i].length = 0;

    //Allocate the particles in the corresponding voxels for space partitioning (neighbors)
    for (let v = 0; v < particles.length; v ++) {
        let particle = particles[v];
        let index = Math.floor(particle.x) + 128 * Math.floor(particle.y);
        if (index >= 0 && index < totalVoxels) voxels[index].push(v);
    }

    let gravity = 10;

    //Predict velocities and positions.
    for(let particle of particles) {
        particle.velocity.y += deltaTime * gravity;
        particle.x += particle.velocity.x * deltaTime;
        particle.y += particle.velocity.y * deltaTime;
    }

    //Solve the iterative position constrain
    for(let i = 0; i < iterations; i ++) iterativeStep();

    //Update the positions
    for(let particle of particles) {
        particle.velocity.x = (particle.x - particle.position.x) / deltaTime;
        particle.velocity.y = (particle.y - particle.position.y) / deltaTime;
        particle.update(particle.x, particle.y);
    }

    //Apply Viscosity
    for(let particle_i of particles) {

        let velocityUpdate = {x:0, y:0}

        for (let u = -1; u <= 1; u++) {
            for (let w = -1; w <= 1; w++) {

                let index = Math.floor(particle_i.x) + u + 128 * (Math.floor(particle_i.y) + w);

                for (let v = 0; v < Math.min(voxels[index].length, maxParticlesPerVoxel); v++) {

                    let particle_j = particles[voxels[index][v]];
                    let distance = particle_i.distance(particle_j);

                    if (distance > 0 && distance < searchRadius) {
                        let partial = searchRadius - distance;
                        velocityUpdate.x += viscosityConstant * (particle_i.velocity.x - particle_j.velocity.x) * partial;
                        velocityUpdate.y += viscosityConstant * (particle_i.velocity.y - particle_j.velocity.y) * partial;
                    }
                }
            }
        }

        particle_i.velocity.x += velocityUpdate.x;
        particle_i.velocity.y += velocityUpdate.y;
    }
}

//simulationStep();

document.body.addEventListener("click", () => {
        looping = !looping;
        if(looping) simulationStep();
    }
);

document.body.addEventListener('keypress', (e) => {
    if(e.key == "n" && !looping) simulationStep();
});
