using Godot;
using System;
using System.Collections.Generic;

public partial class Test : Node2D
{
	// Constants
	private float smoothingLength = 10.0f;
	private float restDensity = 0.59f;
	private float gasConstant = 30.0f;
	private float viscosityCoefficient = 0.1f;
	

	// Particle properties
	private float particleMass = 0.02f;
	private float particleRadius = 3.0f;
	private float particleEffectRadius = 8.0f;
	private List<Particle> particles = new List<Particle>();

	// Simulation parameters
	private Vector2 gravity = new Vector2(0, 9.8f);
	private float timeStep = 0.1f;

	private float cornerTop = 0;
	private float cornerBottom = 40;
	private float cornerLeft = 0;
	private float cornerRight = 300;
	private float spawnXOffset = 7;
	private float spawnYOffset = 7;
	private float spawnXMax = 30;
	private float spawnYMax = 5;
	float M_PI = 3.14159265358979323846f;

	// Initialization
	public override void _Ready()
	{
		GD.Seed(12345);
		RandomNumberGenerator random = new RandomNumberGenerator();
		random.Randomize();
		float volume = particleEffectRadius * particleEffectRadius * M_PI;
		// Create particles
		for (int i = 0; i < spawnXMax; i++)
		{
			for (int j = 0; j < spawnYMax; j++)
			{
				Particle particle = new Particle();
				particle.gravity = gravity;
				particle.Position = new Vector2(i * spawnXOffset+random.RandfRange(0, 5), j*spawnYOffset+random.RandfRange(0, 5));
				particle.Mass = volume * restDensity;
				//particle.Mass = particleMass;
				particle.m_Effect_radius = particleEffectRadius;
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
		DrawLine(new Vector2(cornerLeft,cornerTop), new Vector2(cornerRight,cornerTop), Colors.Red);
		DrawLine(new Vector2(cornerLeft,cornerTop), new Vector2(cornerLeft,cornerBottom), Colors.Red);
		DrawLine(new Vector2(cornerLeft,cornerBottom), new Vector2(cornerRight,cornerBottom), Colors.Red);
		DrawLine(new Vector2(cornerRight,cornerTop), new Vector2(cornerRight,cornerBottom), Colors.Red);
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

		public float				m_Effect_radius = 20f;
		float				m_minRadius;
	
		float				m_stiffness = 50.0f;
		float				m_particleRadiusRatio = 3.0f;
		float				m_viscosity = 0.1f;
		float				m_maxSpeed = 1000.0f;
		float				m_maxAcceleration = 9000.0f;
		float				m_timeScale = 1.0f; // use this to make simulation more stable
		float				m_wallFriction = 0.4f;
		float				m_wallRestitution = 0.8f;
		float               M_PI = 3.14159265358979323846f;
		bool 			    m_bPressure = true;
		bool 			    m_bViscosity = true;

		public bool IsNanOrInfinity(float f)
		{
			if(float.IsNaN(f) || float.IsInfinity(f))
			{
				GD.PrintErr("IsNanOrInfinity Float");
				return true ;
			}
				
			return false;
		}
		public bool IsNanOrInfinity(Vector2 v)
		{
			if(IsNanOrInfinity(v.X) || IsNanOrInfinity(v.Y))
				{
					GD.PrintErr("IsNanOrInfinity Vector2");
					return true ;
				}
			return false;
		}
		
		// Compute density and pressure of the particle
		public void ComputeDensityAndPressure(List<Particle> particles, float smoothingLength, float particleMass, float restDensity, float gasConstant)
		{
			Density = 0.0f;
			float radius = m_Effect_radius;
			float baseWeight = KernelPoly6Density(0.0f, radius);
			

			Density = baseWeight;
				
			
			
			foreach (Particle particle in particles)
			{
				if(particle == this)
				{
					continue;
				}
				float length = Position.DistanceTo(particle.Position);
				float r = Position.DistanceTo(particle.Position);
				if(IsNanOrInfinity(r))
				{
					r=0;
					continue;
				}
				
				
				if(length >= m_Effect_radius)
				{
					continue;
				}
		
				{
					float q = r / m_Effect_radius;
					float w = KernelPoly6Density(r,radius);
					if(IsNearZero(w) || IsNanOrInfinity(w))
						w=0.0001f;
					Density += w;
					
				}
				
				
			}
			Density *= Mass;
			
			Pressure = m_stiffness * (Density - restDensity);
			//GD.Print(Pressure);

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
				
				if(IsNearZero((2.0f * Density * particle.Density)))
					continue;
				if(m_bPressure )
				{
					//Pressure Gradient Forces
					Vector2 pressureAcc = r * (0-Mass) * ((Pressure + particle.Pressure) / (2.0f * Density * particle.Density)) * KernelSpikyGradientFactor(length, m_Effect_radius);
					pressureAcc += r * 0.02f * Mass * ((m_stiffness * (Density + particle.Density)) / (2.0f * Density + particle.Density)) * KernelSpikyGradientFactor(length * 0.8f, m_Effect_radius);
					
					Forces += pressureAcc;
					particle.Forces -= pressureAcc;
				}
				if(m_bViscosity)
				{
				

					{
						
						Vector2 deltaVel = Velocity - particle.Velocity;
						Vector2 viscosityAcc = deltaVel * -Mass * (viscosityCoefficient / (2.0f * Density * particle.Density)) * KernelViscosityLaplacian(length, m_Effect_radius);

						Forces += viscosityAcc;
						particle.Forces -= viscosityAcc;
					}
				}
				
				//particle.Forces -= pressureAcc;


			}
		}
		public bool IsNearZero(float val)
		{
			return (val < 0.0001f && val > -0.0001f);
		}
		// Integrate particle's position and velocity
		public void Integrate(float timeStep)
		{
			Velocity += timeStep * Forces / Density;
			//if(Velocity.Length() > m_maxSpeed)
			{
				//Velocity = Velocity.Normalized() * m_maxSpeed;
			}
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

		float KernelPoly6Density(float r, float h)
		{
			float h2 = h * h;
			float h3 = h2 * h;
			float h9 = h3 * h3* h3;
			float kernel = h2 - r * r;
			float hr =(h2-r*r);
			return (315f*(hr*hr*hr))/(64f*(float)M_PI*h9);
		}
		public void IsOutBox(float cornerLeft, float cornerRight, float cornerTop, float cornerBottom)
		{
			if(IsNanOrInfinity(Position) || IsNanOrInfinity(Velocity))
			{
				Position = new Vector2(cornerLeft,cornerBottom);
				Velocity = new Vector2(0,0);
				GD.PrintErr("ssssssssszzz");
			}
			if (Position.X < cornerLeft)
			{
				Position = new Vector2(cornerLeft,Position.Y);
				Velocity = new Vector2(-Velocity.X*m_wallRestitution,Velocity.Y);
			}
			if (Position.X > cornerRight)
			{
				Position = new Vector2(cornerRight, Position.Y);
				Velocity = new Vector2(-Velocity.X*m_wallRestitution, Velocity.Y);
				//Velocity = new Vector2(0,0);
			}
			if (Position.Y < cornerTop)
			{
				Position =new  Vector2(Position.X, cornerTop);
				Velocity = new Vector2(Velocity.X, -Velocity.Y*m_wallRestitution);
				//Velocity = new Vector2(0,0);
			}
			if (Position.Y > cornerBottom)
			{
				Position = new Vector2(Position.X, cornerBottom);
				Velocity = new Vector2(Velocity.X, -Velocity.Y*m_wallRestitution);
				//Velocity = new Vector2(0,0);
			}
			
		}
		public void AddRandomEpsilonOffsetToPos()
		{
			Position += new Vector2((float)GD.RandRange(-0.0001f, 0.0001f), (float)GD.RandRange(-0.0001f, 0.0001f));
		}

	}

	

	// Main update loop
	public override void _Process(double delta)
	{
		// Compute densities and pressures
		foreach (Particle particle in particles)
		{
			particle.AddRandomEpsilonOffsetToPos();
		}
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
		//Listen For input Key Press Arrow left
		if (Input.IsActionPressed("ui_left"))
		{
			//Move the box to the left
			cornerLeft -= 10;
			cornerRight -= 10;
		}
		//Listen For input Key Press Arrow right
		if (Input.IsActionPressed("ui_right"))
		{
			//Move the box to the right
			cornerLeft += 10;
			cornerRight += 10;
		}
		//Listen For input Key Press Arrow up
		if (Input.IsActionPressed("ui_up"))
		{
			//Move the box to the up
			cornerTop -= 10;
			cornerBottom -= 10;
		}
		//Listen For input Key Press Arrow down
		if (Input.IsActionPressed("ui_down"))
		{
			//Move the box to the down
			cornerTop += 10;
			cornerBottom += 10;
		}

		QueueRedraw();
	}
}
