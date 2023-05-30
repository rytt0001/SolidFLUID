using Godot;
using System;
using System.Collections.Generic;

public partial class Test : Node2D
{
	// Constants
	private float smoothingLength = 10.0f;
	private float restDensity = 1000.0f;
	private float gasConstant = 2000.0f;
	private float viscosityCoefficient = 100.0f;
	

	// Particle properties
	private float particleMass = 1.0f;
	private float particleRadius = 5.0f;
	private List<Particle> particles = new List<Particle>();

	// Simulation parameters
	private Vector2 gravity = new Vector2(0, 9.8f);
	private float timeStep = 0.01f;

	private float cornerTop = 0;
	private float cornerBottom = 600;
	private float cornerLeft = 0;
	private float cornerRight = 800;


	// Initialization
	public override void _Ready()
	{
		// Create particles
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				Particle particle = new Particle();
				particle.gravity = gravity;
				particle.Position = new Vector2(i * 20, j*20);
				particle.Mass = particleMass;
				particles.Add(particle);
			}
			
		}
	}
	// Draw particle
	public override void _Draw()
	{
		foreach (Particle particle in particles)
		{
			DrawParticle(particle);
		}
	}
	public void DrawParticle(Particle p)
	{
		DrawCircle(p.Position, particleRadius, Colors.White);
		DrawLine(new Vector2(cornerLeft,cornerTop), new Vector2(cornerRight,cornerTop), Colors.Red,2);
		DrawLine(new Vector2(cornerLeft,cornerTop), new Vector2(cornerLeft,cornerBottom), Colors.Red,2);
		DrawLine(new Vector2(cornerLeft,cornerBottom), new Vector2(cornerRight,cornerBottom), Colors.Red,2);
		DrawLine(new Vector2(cornerRight,cornerTop), new Vector2(cornerRight,cornerBottom), Colors.Red,2);
	}
	// Particle class
	public class Particle
	{
		public Vector2 Position { get; set; }
		public Vector2 Velocity { get; set; }
		public float Density { get; set; }
		public float Pressure { get; set; }
		public float Mass { get; set; }
		public Vector2 Forces { get; set; }
		public Vector2 gravity { get; set; }

		float				m_Effect_radius = 20f;
		float				m_minRadius;
		float				m_restDensity = 0.59f;
		float				m_stiffness = 500.0f;
		float				m_particleRadiusRatio = 3.0f;
		float				m_viscosity = 0.1f;
		float				m_maxSpeed = 10.0f;
		float				m_maxAcceleration = 900.0f;
		float				m_timeScale = 1.0f; // use this to make simulation more stable
		float				m_wallFriction = 0.4f;
		float				m_wallRestitution = 0.4f;
		float               M_PI = 3.14159265358979323846f;
		// Compute density and pressure of the particle
		public void ComputeDensityAndPressure(List<Particle> particles, float smoothingLength, float particleMass, float restDensity, float gasConstant)
		{
			Density = 0.0f;
			foreach (Particle particle in particles)
			{
				float r = Position.DistanceTo(particle.Position);
				if (r < smoothingLength)
				{
					float q = r / smoothingLength;
					float w = KernelCubicSpline(q);
					Density += particleMass * w;
				}
			}

			Pressure = gasConstant * (Density - restDensity);
		}

		// Compute forces acting on the particle
		public void ComputeForces(List<Particle> particles, float smoothingLength, float particleMass, float viscosityCoefficient)
		{
			Forces = gravity * Density;

			foreach (Particle particle in particles)
			{
				if(particle == this)
				{
					continue;
				}
				float length = Position.DistanceTo(particle.Position);
				if(length >= m_Effect_radius)
				{
					continue;
				}
				Vector2 r = Position - particle.Position;
				


				//Pressure Gradient Forces
				Vector2 pressureAcc = r * (0-Mass) * ((Pressure + particle.Pressure) / (2.0f * Density * particle.Density)) * KernelSpikyGradientFactor(length, m_Effect_radius);
				pressureAcc += r * 0.02f * Mass * ((m_stiffness * (Density + particle.Density)) / (2.0f * Density + particle.Density)) * KernelSpikyGradientFactor(length * 0.8f, m_Effect_radius);
				Forces += pressureAcc;
				GD.Print(KernelSpikyGradientFactor(length, m_Effect_radius));
				//particle.Forces -= pressureAcc;


			}
		}

		// Integrate particle's position and velocity
		public void Integrate(float timeStep)
		{
			Velocity += timeStep * Forces / Density;
			Position += timeStep * Velocity;
		}

		

		// Cubic spline kernel function
		private float KernelCubicSpline(float q)
		{
			float result = 0.0f;
			if (q >= 0.0f && q <= 1.0f)
			{
				float factor = 2.0f / 3.0f;
				float term1 = 2.0f / 3.0f - q * q + 0.5f * q * q * q;
				float term2 = 1.0f - 1.5f * q * q + 0.75f * q * q * q;
				result = factor * (term1 * term1 * term1 - term2 * term2 * term2);
			}

			return result;
		}

		float KernelDefault(float r, float h)
		{
			float h2 = h * h;
			float h4 = h2 * h2;
			float kernel = h2 - r * r;
			return (kernel * kernel * kernel) * (4.0f / (((float)M_PI) * h4 * h4));
		}

		float KernelDefaultGradientFactor(float r, float h)
		{
			float h2 = h * h;
			float h4 = h2 * h2;
			float kernel = h2 - r * r;
			return -(kernel * kernel) * (6.0f / (((float)M_PI) * h4 * h4));
		}

		float KernelDefaultLaplacian(float r, float h)
		{
			float h2 = h * h;
			float h4 = h2 * h2;
			float kernel = h2 - r * r;
			return -(kernel * kernel) * (6.0f / (((float)M_PI) * h4 * h4)) * (3.0f * h2 - 7.0f * r * r);
		}


		float KernelSpikyGradientFactorNorm(float r, float h)
		{
			float h2 = h * h;
			float h5 = h2 * h2 * h;
			float kernel = h - r;
			return kernel * kernel * (-15.0f / ((float)M_PI * h5));
		}

		float KernelSpikyGradientFactor(float r, float h)
		{
			float h2 = h * h;
			float h5 = h2 * h2 * h;
			float kernel = h - r;
			
			return kernel * kernel * (-15.0f / ((float)M_PI * h5 * r));
		}

		float KernelViscosityLaplacian(float r, float h)
		{
			float h2 = h * h;
			float kernel = h - r;
			return kernel * (30.0f / ((float)M_PI * h2 * h2 * h));
		}

		float KernelPoly6hGradientFactor(float r, float h)
		{
			float h2 = h * h;
			float kernel = h2 - r * r;
			return kernel * kernel * (24.0f / ((float)M_PI * h2 * h2 * h2 * h * r));
		}

		// Compute gradient of the SPH spiky kernel
		private Vector2 GradientSpikyKernel(Vector2 r, float smoothingLength)
		{
			float rLength = r.Length();
			if (rLength > 0 && rLength < smoothingLength)
			{
				float factor = -45.0f / (3.14f * Mathf.Pow(smoothingLength, 6));
				float scalar = smoothingLength - rLength;
				return factor * scalar * scalar * r.Normalized();
			}
			else
			{
				return Vector2.Zero;
			}
		}

		public void IsOutBox(float cornerLeft, float cornerRight, float cornerTop, float cornerBottom)
		{
			if (Position.X < cornerLeft)
			{
				Position = new Vector2(cornerLeft,Position.Y);
				Velocity = new Vector2(-Velocity.X,Velocity.Y);
			}
			if (Position.X > cornerRight)
			{
				Position = new Vector2(cornerRight, Position.Y);
				Velocity = new Vector2(-Velocity.X, Velocity.Y);
				Velocity = new Vector2(0,0);
			}
			if (Position.Y < cornerTop)
			{
				Position =new  Vector2(Position.X, cornerTop);
				Velocity = new Vector2(Velocity.X, -Velocity.Y);
				Velocity = new Vector2(0,0);
			}
			if (Position.Y > cornerBottom)
			{
				Position = new Vector2(Position.X, cornerBottom);
				Velocity = new Vector2(Velocity.X, -Velocity.Y);
				Velocity = new Vector2(0,0);
			}
			
		}

	}

	// Main update loop
	public override void _Process(double delta)
	{
		// Compute densities and pressures
		foreach (Particle particle in particles)
		{
			particle.ComputeDensityAndPressure(particles, smoothingLength, particleMass, restDensity, gasConstant);
		}

		// Compute forces
		foreach (Particle particle in particles)
		{
			particle.ComputeForces(particles, smoothingLength, particleMass, viscosityCoefficient);
		}

		// Integrate positions and velocities
		foreach (Particle particle in particles)
		{
			particle.IsOutBox(cornerLeft, cornerRight, cornerTop, cornerBottom);
			particle.Integrate(timeStep);
		}
		

		QueueRedraw();
	}
}
