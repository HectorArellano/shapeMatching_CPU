
/*
Shape Matching for Position Based Dynamics
based on the paper from Matthias Mueller:
 http://matthias-mueller-fischer.ch/talks/2017-EG-CourseNotes.pdf
 */


import {Particle} from './Particle.js';

const radius = 6; //Visual radius of the particles
const deltaTime = 0.06; //speed of the simulation
const center = {x:64, y:64}; //center of the border circle
const circleRadius = 64 * 0.494; //radius of the circle
const iterations = 30; //Number of iterations for the shape matching solver
const gravity = 10;


let shapes = [];
let particles = [];
let ii = 0;
let container = document.querySelector('.container');
let looping = true;

let voxels = [];
let totalVoxels = 128 * 128;
let globalParticles =[];
for (let i = 0; i < totalVoxels; i++) voxels[i] = [];



function evaluateCenterOfMass(particles) {
    let centerOfMass = {x: 0, y:0}
    for(let particle of particles) {
        centerOfMass.x += particle.position.x;
        centerOfMass.y += particle.position.y;
    }
    centerOfMass.x /= particles.length;
    centerOfMass.y /= particles.length;
    return centerOfMass;
}

let upperBorder = document.querySelector(".border2");
upperBorder.x = 0;
upperBorder.style.left = String(upperBorder.x).concat("px");

document.addEventListener("mousemove", (e) => {
    upperBorder.x = e.clientX;
    upperBorder.style.left = String(upperBorder.x).concat("px");
})

/*
The name in the shape and the particles is used to allocate the particles in the DOM (to place them in their positions),
but most important it's used for collisions, where a single pool of particles will be used to define the neighborhood
search for close collisions, the particles from the same shape shouldn't be evaluated for collisions (avoiding inner
collisions).

cx... center X
cy... center Y
stiffness... controls the soft body softness, the range goes from [0 - 1], values above 0.4 are advised.
 */

let shapeDefinition = [
    {radius: 7, cx: 50, cy: 50, stiffness: 0.3},
    {radius: 7, cx: 65, cy: 50, stiffness: 0.3},
    {radius: 7, cx: 80, cy: 50, stiffness: 0.3}
]


//create the particles and shapes.
for(let shapeDef of shapeDefinition) {
    particles = [];

    //Generate the voxelized particles allocated into the shape / volume
    //For a general mesh this should be generated using a voxelization process (meshes should be closed).
    for(let j = 0; j < 128; j ++) {
        for (let i = 0; i < 128; i++) {
            let r = Math.pow((i - shapeDef.cx), 2) + Math.pow((j - shapeDef.cy), 2);
            if(r < shapeDef.radius * shapeDef.radius) {
                let p = new Particle(`shape${ii}`, container, i, j, particles.length);
                particles.push(p);
                particles[particles.length - 1].radius(radius);
                globalParticles.push(p);
            }
        }
    }

    //actual shape, can be a regular mesh.
    for(let i = 0; i < 60; i ++) {
        const angle = 2 * Math.PI * i / 60;
        const r = shapeDef.radius;
        let p = new Particle(`shape${ii}`, container, r * Math.cos(angle) + shapeDef.cx, r * Math.sin(angle) + shapeDef.cy, particles.length);
        particles.push(p);
        particles[particles.length - 1].radius(radius);
        globalParticles.push(p);
    }

    shapes[ii] = {radius: Math.sqrt(90), name: `shape${ii}`, particles: [...particles], initialCenterOfMass: evaluateCenterOfMass(particles), stiffness: shapeDef.stiffness, centerOfMass: evaluateCenterOfMass(particles)}

    ii++;
}

//TODO this is precalculated should be saved in a texture. Should check if it's better to generate this in the CPU or the GPU
//Evaluate the relative position for each particle relative to the center of mass to each shape
for(let shape of shapes) {

    //Evaluate Aqq for the linear extension
    let Aqq = mat2.fromValues(0, 0, 0, 0);

    for(let particle of shape.particles) {
        particle.initialRelativePosition.x = particle.position.x - shape.initialCenterOfMass.x;
        particle.initialRelativePosition.y = particle.position.y - shape.initialCenterOfMass.y;

        const q = particle.initialRelativePosition;
        let pqi = mat2.fromValues(q.x * q.x, q.y * q.x, q.x * q.y, q.y * q.y);

        mat2.add(Aqq, Aqq, pqi);
    }

    mat2.invert(Aqq, Aqq);
    shape.Aqq = Aqq;
}

particles = null;



/*
 Collision evaluation between particles and borders
 This is done particle wise and it could be O(n2) without
 any acceleration:
 - Acceleration is done evaluating bounding volumes / shapes
 - For the GPU a neighborhood search could be implemented  <<--- implies Harada sorting in buckets.
 */

let evaluateCollisions = (shape) => {

    for(const particle of shape.particles) {

        const r = 0.7;

        //Grid search
        for (let u = -2; u <= 2; u++) {
            for (let w = -2; w <= 2; w++) {

                let index = Math.floor(particle.x) + u + 128 * (Math.floor(particle.y) + w);

                for(let v = 0; v < Math.min(voxels[index].length, 100); v ++) {

                    let c_particle = globalParticles[voxels[index][v]];

                    if(c_particle.name != particle.name) {
                        let normal = {x: particle.position.x - c_particle.position.x, y: particle.position.y - c_particle.position.y}
                        let normalL = Math.sqrt(normal.x * normal.x + normal.y * normal.y);

                        //Move the particle to the outer radius of the two that are in contact.
                        if (normalL < 2 * r && normalL > 0) {

                            normal.x /= normalL;
                            normal.y /= normalL;

                            particle.position.x += r * normal.x;
                            particle.position.y += r * normal.y;

                            //TODO check how to deal with the symmetry in the GPU
                            c_particle.position.x -= r * normal.x;
                            c_particle.position.y -= r * normal.y;
                        }
                    }
                }
            }
        }

//        //Collision among particles O(n2)
//        for (const c_shape of shapes) {
//
//            if (c_shape.name != particle.name) {
//
//                let normal = {x: particle.position.x - c_shape.centerOfMass.x, y: particle.position.y - c_shape.centerOfMass.y}
//                let normalL = Math.sqrt(normal.x * normal.x + normal.y * normal.y);
//
//                //There's a collision with the bounding shape (circle), should evaluate the inner particles.
//                if(normalL < 1.5 * c_shape.radius) {
//
//                    //Check the inner particles of the shape
//                    for (const c_particle of c_shape.particles) {
//
//                        normal = {x: particle.position.x - c_particle.position.x, y: particle.position.y - c_particle.position.y}
//                        normalL = Math.sqrt(normal.x * normal.x + normal.y * normal.y);
//
//                        //Move the particle to the outer radius of the two that are in contact.
//                        if (normalL < 2 * r && normalL > 0) {
//
//                            normal.x /= normalL;
//                            normal.y /= normalL;
//
//                            particle.position.x += r * normal.x;
//                            particle.position.y += r * normal.y;
//
//                            //TODO check how to deal with the symmetry in the GPU
//                            c_particle.position.x -= r * normal.x;
//                            c_particle.position.y -= r * normal.y;
//                        }
//                    }
//                }
//            }
//        }

        let normal = {x: particle.position.x - center.x, y: particle.position.y - center.y}
        let normalL = Math.sqrt(normal.x * normal.x + normal.y * normal.y);
        var dist = normalL - circleRadius;
        var lMax = 0.0001;
        if (dist > 0.) {
            normal.x /= normalL;
            normal.y /= normalL;
            var cP = {x: center.x + circleRadius * normal.x, y: center.y + circleRadius * normal.y}
            particle.position.x = cP.x - lMax * normal.x;
            particle.position.y = cP.y - lMax * normal.y;
            particle.prevPosition.x = particle.position.x;
            particle.prevPosition.y = particle.position.y;
        }


        if (particle.position.x < (upperBorder.x + 20) / 8) {
            particle.position.x = (upperBorder.x + 20) / 8;
            particle.prevPosition.x = particle.position.x;
        }
    }
}



/*
This function is used to evaluate the rotation matrix used to update the particles
it's called per shape, but requires to sum all the particles.
 */
let generateShapeMatrix = shape => {

    //Evaluate the relative position for each particle relative to the center of mass O(n)
    let Apq = mat2.fromValues(0, 0, 0, 0);
    for(let particle of shape.particles) {
        particle.relativePosition.x = particle.position.x - shape.centerOfMass.x;
        particle.relativePosition.y = particle.position.y - shape.centerOfMass.y;

        const p = particle.relativePosition;
        const q = particle.initialRelativePosition;
        let pqi = mat2.fromValues(p.x * q.x, p.y * q.x, p.x * q.y, p.y * q.y);
        mat2.add(Apq, Apq, pqi);
    }

    //Calculate the transpose of Apq
    let ApqT = mat2.create();
    mat2.transpose(ApqT, Apq);

    //Calculate S
    let S = mat2.create();
    mat2.mul(S, ApqT, Apq);

    /*
     Using Babylonian method to evaluate the square root.
     */

    let Si = mat2.create();
    let SiI = mat2.create();
    for(let i = 0; i < 50; i ++) {
        SiI = mat2.invert(SiI, Si);
        mat2.multiplyScalar(Si, mat2.add(Si, Si, mat2.multiply(SiI, S, SiI)), 0.5);
    }

    let Sinv = mat2.create();
    mat2.invert(Sinv, Si);

    //Calculate rotation matrix
    let R = mat2.create();
    mat2.mul(R, Apq, Sinv);


    //Linear extension
    let A = mat2.create();
    A = mat2.mul(A, Apq, shape.Aqq);
    let detA = mat2.determinant(A);
    detA = Math.pow(detA, 1/ 3);
    mat2.multiplyScalar(A, A, detA);

    const beta = 0.5;
    mat2.add(A, mat2.multiplyScalar(A, A, beta), mat2.multiplyScalar(R, R, 1 - beta));

    return A;
}



//Calculate the displacement for all the particles of the shape
let updateParticlesIterationPositions = (shape, A) => {
    for (let particle of shape.particles) {
        const d = particle.initialRelativePosition;
        let g = vec2.fromValues(d.x, d.y);

        //Apply the rotation for the initial relative position.
        vec2.transformMat2(g, g, A);

        particle.position.x += shape.stiffness * (g[0] + shape.centerOfMass.x - particle.position.x);
        particle.position.y += shape.stiffness * (g[1] + shape.centerOfMass.y - particle.position.y);

    }
}



let iterativeStep = (shape) => {

    //TODO GPU this should be a draw call on a small texture that updates the data for the shapes CM
    for(let shape of shapes) shape.centerOfMass = evaluateCenterOfMass(shape.particles); //<<<<<<<--------- O(n)


    //TODO GPU possible bottleneck rendering all the particles from the shape.
    evaluateCollisions(shape); //<<<<---------- O(n2)


    //TODO GPU draw call to update the A matrix used to rotate
    const A = generateShapeMatrix(shape); //<<<<<<<<<--------- O(n)


    //TODO GPU draw call to update the iteration position of each particle
    updateParticlesIterationPositions(shape, A);
}

let simulationStep = () => {

    if(looping) requestAnimationFrame(simulationStep);


    //Set the particles inside the voxels
    for (let i = 0; i < totalVoxels; i++) voxels[i].length = 0;

    //Allocate the particles in the corresponding voxels for space partitioning (neighbors)
    for (let v = 0; v < globalParticles.length; v ++) {
        let particle = globalParticles[v];
        let index = Math.floor(particle.x) + 128 * Math.floor(particle.y);
        if (index >= 0 && index < totalVoxels) voxels[index].push(v);
    }


    for (let shape of shapes) {

        //TODO GPU write a shader that updates the iteration positions for the step
        for(let particle of shape.particles) {
            particle.velocity.y += deltaTime * gravity;
            particle.position.x += particle.velocity.x * deltaTime;
            particle.position.y += particle.velocity.y * deltaTime;
        }


        //Solve the iterative position constrain
        //TODO GPU write the corresponding draw calls for the iterative process.
        for(let i = 0; i < iterations; i ++) iterativeStep(shape);



        //TODO GPU write a shader that updates the final positions for the step
        for(let particle of shape.particles) {

            particle.x = particle.position.x;
            particle.y = particle.position.y;

            particle.velocity.x = (particle.x - particle.prevPosition.x) / deltaTime;
            particle.velocity.y = (particle.y - particle.prevPosition.y) / deltaTime;

            particle.update(particle.x, particle.y);
        }
    }
}

simulationStep();

document.body.addEventListener("click", () => {
        looping = !looping;
        if(looping) simulationStep();
    }
);

document.body.addEventListener('keypress', (e) => {
    if(e.key == "n" && !looping) simulationStep();
});
