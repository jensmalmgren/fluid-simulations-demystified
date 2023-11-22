using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace FluidSimulationsDemystified
{
	public partial class Form : System.Windows.Forms.Form
	{
		// https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf
		// https://youtu.be/qsYE1wMEMPA?si=QhRSunPvcCLazm70

		public DrawPanel _DrawPanel;
		private Fluid f;
		Random m_Random = new Random(137);
		public Form()
		{
			InitializeComponent();

			Width = (fa.NX + 2) * fa.displayScale;
			Height = (fa.NY + 4) * fa.displayScale;
			Text = "Fluid simulations by Jos Stam demystified by Jens Malmgren";

			_DrawPanel = new DrawPanel();
			_DrawPanel.Dock = DockStyle.Fill;
			Controls.Add(_DrawPanel);

			f = new Fluid();

			_DrawPanel.Paint += OnDrawPanel_Paint;
		}

		private void OnDrawPanel_Paint(object sender, PaintEventArgs e)
		{
			float _fDensity = 1f;
			float _fVelocity = 400f;

			fa.angle += (float)m_Random.Next(-10, 11);

			if (fa.angle > 360f)
			{
				fa.angle = fa.angle - 360;
			}
			else
			{
				if (fa.angle < 0f)
				{
					fa.angle = fa.angle + 360;
				}
			}

			float _xDir = _fVelocity * (float)Math.Sin(fa.angle * Math.PI / 180f);
			float _yDir = _fVelocity * (float)Math.Cos(fa.angle * Math.PI / 180f);
			for (int i = 0; i < 1; i++)
			{
				for (int j = 0; j < 1; j++)
				{
					f.d1[i + fa.NX / 2, j + fa.NY / 2] = _fDensity;
					f.set_velocity(i + fa.NX / 2, j + fa.NY / 2, _xDir, _yDir);
				}
			}

			for (int x = 0; x < fa.NX; x++)
			{
				for (int y = 0; y < fa.NY; y++)
				{
					float _density = f.d1[x, y];
					float _ratioWithMaxDensity = _density / _fDensity;
					int _c = fa.Clamp((int)(_ratioWithMaxDensity * 255f) , 0, 255);
					Color _Color = Color.FromArgb(_c, _c, _c);
					using (SolidBrush _Solid = new SolidBrush(_Color))
					{
						e.Graphics.FillRectangle(_Solid, x * fa.displayScale, y * fa.displayScale, fa.displayScale, fa.displayScale);
					}
				}
			}

			f.sim_step();
			_DrawPanel.Invalidate();
		}
	} // class Form

	/// <summary>
	/// "Fluid All" global variables.
	/// </summary>
	public static class fa
	{
		public static int NX = 50;
		public static int NY = 50;
		public static int iter = 4;
		public static int displayScale = 9;
		public static float angle = 0;
		public static float deltat = 0.01f;
		public static float diff = 0.00001f;
		public static int iTraceColumnWidth = 6;
		public static bool bSetBoundary = false;

		public static int Clamp(int value, int min, int max)
		{
			return (value < min) ? min : (value > max) ? max : value;
		}
	}

	public enum Boundary
	{
		None = 0,
		LeftAndRight = 1,
		TopAndBottom = 2
	}

	public class Fluid
	{
		public Matrix u1 = new Matrix("u1");
		public Matrix u2 = new Matrix("u2");
		public Matrix u3 = new Matrix("u3");
		public Matrix u4 = new Matrix("u4");
		public Matrix u4tmp = new Matrix("u4tmp");
		public Matrix u5 = new Matrix("u5");
		public Matrix u6 = new Matrix("u6");
		public Matrix u7 = new Matrix("u7");
		public Matrix u8 = new Matrix("u8");
		public Matrix u9 = new Matrix("u9");
		public Matrix u10 = new Matrix("u10");
		public Matrix v1 = new Matrix("v1");
		public Matrix v2 = new Matrix("v2");
		public Matrix v3 = new Matrix("v3");
		public Matrix v4 = new Matrix("v4");
		public Matrix v4tmp = new Matrix("v4tmp");
		public Matrix v5 = new Matrix("v5");
		public Matrix v6 = new Matrix("v6");
		public Matrix v7 = new Matrix("v7");
		public Matrix v8 = new Matrix("v8");
		public Matrix v9 = new Matrix("v9");
		public Matrix v10 = new Matrix("v10");
		public Matrix p1 = new Matrix("p1");
		public Matrix p2 = new Matrix("p2");
		public Matrix p3 = new Matrix("p3");
		public Matrix p4 = new Matrix("p4");
		public Matrix p5 = new Matrix("p5");
		public Matrix p6 = new Matrix("p6");
		public Matrix p7 = new Matrix("p7");
		public Matrix p8 = new Matrix("p8");
		public Matrix p9 = new Matrix("p9");
		public Matrix p10 = new Matrix("p10");
		public Matrix div1 = new Matrix("div1");
		public Matrix div2 = new Matrix("div2");
		public Matrix d1 = new Matrix("d1");
		public Matrix d2 = new Matrix("d2");
		public Matrix d3 = new Matrix("d3");
		public Matrix d4 = new Matrix("d4");
		public Matrix d4tmp = new Matrix("d4tmp");
		public Matrix d5 = new Matrix("d5");
		public Matrix d6 = new Matrix("d6");
		public Matrix d7 = new Matrix("d7");
		public Matrix d8 = new Matrix("d8");
		public Matrix d9 = new Matrix("d9");

		public void set_velocity(int x, int y, float amountX, float amountY)
		{
			u1[x, y] = fa.deltat * amountX;
			v1[x, y] = fa.deltat * amountY;
		} // set_velocity()

		public void sim_step()
		{
			// In Jos Stam's code, this part was called vel_step.

			// Jos Stam used the suffix 0 to indicate that a matrix
			// held the previous values. We need numbers for other things
			// so that is inconvenient. Besides, his program is
			// rewritten by me in such a way no matrix is reused.
			// In that scenario, current and previous is meaningless.

			// u stands for horizontal velocity vectors and v for vertical velocity vectors.

			// Here, we calculate based on vectors 1 and 2 and leave the result in 3.
			add_source(u1, u2, u3);
			add_source(v1, v2, v3);

			// In the original paper by Jos Stam, u and u0 were swapped
			// with as result that the values of u went to u0.
			// We swap by argument order, creating the same result
			// with a copy. Stams u after swap is our u4. Stams u0 is our u3.

			u3.Copy(u4);
			u3.Copy(u4tmp);

			// Jos Stam is saying that diffusion is calculated with the
			// mean value in a star-shaped form of adjacent cells.
			// That is true, but since the calculation is moving
			// forward in steps of one, the next star is overlapping
			// the previous. Hence, a moving average is created,
			// and I don't think that is intended. Because of this
			// I introduced a temp matrix u4tmp to resolve this.

			// Jos Stam used an iteration of 20 in the diffuse function,
			// but four repeats are sufficient. Since we do all the calculations
			// by handing over the computed result into new matrixes,
			// the iteration is done with four consecutive calls,
			// where u4->u5, u5->u6 and so on.

			diffuse_noiter(Boundary.LeftAndRight, u4, u3, u4tmp, u5);
			diffuse_noiter(Boundary.LeftAndRight, u5, u3, u4tmp, u6);
			diffuse_noiter(Boundary.LeftAndRight, u6, u3, u4tmp, u7);
			diffuse_noiter(Boundary.LeftAndRight, u7, u3, u4tmp, p1);

			// We prepare for the function project in the final
			// diffusion step. u7->p1.

			// Here, we are doing the same calculation for the vertical component.

			v3.Copy(v4);
			v3.Copy(v4tmp);

			diffuse_noiter(Boundary.TopAndBottom, v4, v3, v4tmp, v5);
			diffuse_noiter(Boundary.TopAndBottom, v5, v3, v4tmp, v6);
			diffuse_noiter(Boundary.TopAndBottom, v6, v3, v4tmp, v7);
			diffuse_noiter(Boundary.TopAndBottom, v7, v3, v4tmp, v8);

			// For v8, it is necessary to make a copy to div1 since div1 is
			// a target matrix and v8 is a source matrix.

			v8.Copy(div1);

			projectv2(u5, v8, p1, div1, p6, div2, p2, p3, p4, p5);

			// In Jos Stam's code, he swaps u with u0, meaning that u3 is swapped with p1.
			// Also v is swapped with v0 meaning it is swapped with div1 (a copy of v8)

			// advect(u,u0,u0,v0)->advect(d,d0,u,v)->d=u. But for us u is p1,	u0 is p6,	u0 is p6, v0 is div2. Result goes to u10.
			// advect(v,v0,u0,v0)->advect(d,d0,u,v)->d=v. But for us v is div1, v0 is div2, u0 is p6, v0 is div2. Result goes to v10.

			advect(Boundary.LeftAndRight, p1,	p6,		p6,		div2, u9);
			advect(Boundary.TopAndBottom, div1, div2,	p6,		div2, v9);

			projectv2(u9, v9, p6, div2, u10, v10, p7, p8, p9, p10);

			// originally this part was dens_step
			add_source(d1, d2, d3);

			d1.Copy(d4);

			diffuse_noiter(0, d4, d3, d4tmp, d5);
			diffuse_noiter(0, d5, d3, d4tmp, d6);
			diffuse_noiter(0, d6, d3, d4tmp, d7);
			diffuse_noiter(0, d7, d3, d4tmp, d8);

			advect(0, d4, d8, u10, v10, d9);

			// Here, we prepare for the next round by copying all end matrices to the start matrices.
			d9.Copy(d1);
			u10.Copy(u1);
			v10.Copy(v1);
			p6.Copy(u2);
			div2.Copy(v2);
		} // sim_step()

		/// <summary>
		/// This function is used to diffuse both velocity as well as density.
		/// Here is the Gauss-Seidel method used. Called four times.
		/// </summary>
		/// <param name="b"></param>
		/// <param name="xs">Current density or velocity. Source matrix.</param>
		/// <param name="xps">Previous density or velocity. Source matrix.</param>
		/// <param name="xtmps">Temporary x matrix. Source matrix.</param>
		/// <param name="xt">The diffused x. Target matrix.</param>
		public void diffuse_noiter(Boundary b, Matrix xs, Matrix xps, Matrix xtmps, Matrix xt)
		{
			xs.Copy(xt);
			float k = fa.deltat * fa.diff * fa.NX * fa.NY;
			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					xt[i, j] = (xps[i, j] + k * (xtmps[i - 1, j    ] + xtmps[i + 1, j    ]
											   + xtmps[i,     j - 1] + xtmps[i,     j + 1])) / (1 + 4 * k);
				}
			}
			set_bnd(b, xt);
		} // diffuse_noiter()

		public void advect(Boundary b, Matrix ds, Matrix DS, Matrix us, Matrix vs, Matrix dt)
		{

			ds.Copy(dt);

			int i0, j0, i1, j1;
			float x, y;

			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					if (us[i, j] == 0f || vs[i, j] == 0f) continue;

					// Find where the density will come from.
					// Get position of square.
					// Subtract velocity multiplied with dt.
					x = (float)i - fa.deltat * fa.NX * us[i, j];
					y = (float)j - fa.deltat * fa.NY * vs[i, j];

					if (x < 0.5f) x = 0.5f; // Keep x inside.
					if (x > (float)fa.NX + 0.5f) x = (float)fa.NX + 0.5f; // Keep x inside.
					if (y < 0.5f) y = 0.5f; // Keep y inside.
					if (y > (float)fa.NY + 0.5f) y = (float)fa.NY + 0.5f; // Keep y inside.

					i0 = (int)x; // floor part of x. Right cells.
					i1 = i0 + 1; // Left cells.

					j0 = (int)y; // floor part of y. Top cells coordinate.
					j1 = j0 + 1; // Bottom cells coordinate.

					//Linear interpolation by Jos Stam.
					//s1 = x - i0;
					//s0 = 1.0f - s1;
					//t1 = y - j0;
					//t0 = 1.0f - t1;
					//dt[i, j] = s0 * (t0 * DS[i0, j0] + t1 * DS[i0, j1]) + s1 * (t0 * DS[i1, j0] + t1 * DS[i1, j1]);

					float z1 = linear_interpolate(DS[i0, j0], DS[i1, j0], x - i0); // k is fraction part of x
					float z2 = linear_interpolate(DS[i0, j1], DS[i1, j1], x - i0); // k is fraction part of x
					dt[i, j] = linear_interpolate(z1, z2, y - j0); // k is fraction part of y
				}
			}
			set_bnd(b, dt);
		} // advect()

		public float linear_interpolate(float a, float b, float k)
		{
			return a + k * (b - a);
		} // linear_interpolate()

		public void add_source(Matrix xs, Matrix s, Matrix xt)
		{
			for (int i = 0; i < fa.NX; i++)
			{
				for (int j = 0; j < fa.NY; j++)
				{
					xt[i,j] = xs[i,j] + fa.deltat * s[i,j];
				}
			}
		} // add_source()

		/// <summary>
		/// Clearing divergence. Helmholtz decomposition.
		/// Helmholtz Theorem a.k.a. the fundamental theorem of vector calculus.
		/// </summary>
		/// <param name="us">Source matrix for u, the horizontal component.</param>
		/// <param name="vs">Source matrix for v, the vertical component.</param>
		/// <param name="pt">The resulting p matrix. Target matrix.</param>
		/// <param name="divt">Div target matrix.</param>
		/// <param name="ut">Horizontal target matrix.</param>
		/// <param name="vt">Vertical target matrix.</param>
		/// <param name="pAt">Values from pt goes here.</param>
		/// <param name="pBt">Values from pAt goes here.</param>
		/// <param name="pCt">Values from pBt goes here.</param>
		/// <param name="pDt">Valued from pCt goes here. This is the result p matrix.</param>
		public void projectv2(Matrix us, Matrix vs, Matrix pt, Matrix divt, Matrix ut, Matrix vt,
			Matrix pAt, Matrix pBt, Matrix pCt, Matrix pDt)
		{
			projectA(us, vs, divt);
			projectB(pt);
			projectC(divt, pt, pAt);
			projectC(divt, pAt, pBt);
			projectC(divt, pBt, pCt);
			projectC(divt, pCt, pDt);
			projectD(pDt, us, vs, ut, vt);
		} // projectv2()

		/// <summary>
		/// Calculate divergence.
		/// </summary>
		/// <param name="us"></param>
		/// <param name="vs"></param>
		/// <param name="divt"></param>
		public void projectA(Matrix us, Matrix vs, Matrix divt)
		{
			float hx = 1.0f / Math.Max(fa.NX, fa.NY);

			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					divt[i, j] = -0.5f * hx * (us[i + 1, j    ] - us[i - 1, j    ] +
											   vs[i,	 j + 1] - vs[i,		j - 1]);
				}
			}
			set_bnd(0, divt);
		} // projectA()

		/// <summary>
		/// Prepare the p matrix.
		/// </summary>
		/// <param name="pt"></param>
		public void projectB(Matrix pt)
		{
			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					pt[i, j] = 0;
				}
			}
			set_bnd(0, pt);
		} // projectB()

		/// <summary>
		/// Gauss seidel method to solve p. Called four times.
		/// </summary>
		/// <param name="divs"></param>
		/// <param name="ps"></param>
		/// <param name="pt"></param>
		public void projectC(Matrix divs, Matrix ps, Matrix pt)
		{
			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					pt[i, j] = (divs[i, j] + ps[i - 1, j    ] + ps[i + 1, j    ] +
							 			     ps[i,     j - 1] + ps[i,     j + 1]) / 4;
				}
			}
			set_bnd(0, pt);
		} // projectC()

		/// <summary>
		/// Subtract divergence to calculate the curl.
		/// </summary>
		/// <param name="ps"></param>
		/// <param name="us"></param>
		/// <param name="vs"></param>
		/// <param name="ut"></param>
		/// <param name="vt"></param>
		public void projectD(Matrix ps, Matrix us, Matrix vs, Matrix ut, Matrix vt)
		{
			float hx = 1.0f / fa.NX;
			float hy = 1.0f / fa.NY;

			if (fa.bSetBoundary)
			{
				us.Copy(ut);
				vs.Copy(vt);
			}
			for (int i = 1; i <= fa.NX; i++)
			{
				for (int j = 1; j <= fa.NY; j++)
				{
					ut[i, j] = us[i, j] - 0.5f * (ps[i + 1, j    ] - ps[i - 1, j    ]) / hx;
					vt[i, j] = vs[i, j] - 0.5f * (ps[i,     j + 1] - ps[i,     j - 1]) / hy;
				}
			}
			set_bnd(Boundary.LeftAndRight, ut);
			set_bnd(Boundary.TopAndBottom, vt);
		} // projectD()

		void set_bnd(Boundary b, Matrix x)
		{
			if (fa.bSetBoundary)
			{
				for (int i = 1; i < fa.NY - 1; i++)
				{
					x[0, i] = b == Boundary.LeftAndRight ? -x[1, i] : x[1, i];
					x[fa.NX + 1, i] = b == Boundary.LeftAndRight ? -x[fa.NX, i] : x[fa.NX, i];
				}

				for (int i = 1; i < fa.NX - 1; i++)
				{
					x[i, 0] = b == Boundary.TopAndBottom ? -x[i, 1] : x[i, 1];
					x[i, fa.NY + 1] = b == Boundary.TopAndBottom ? -x[i, fa.NY] : x[i, fa.NY];
				}

				x[0, 0] = 0.5f * (x[1, 0] + x[0, 1]);
				x[0, fa.NY + 1] = 0.5f * (x[1, fa.NY + 1] + x[0, fa.NY]);
				x[fa.NX + 1, 0] = 0.5f * (x[fa.NX, 0] + x[fa.NX + 1, 1]);
				x[fa.NX + 1, fa.NY + 1] = 0.5f * (x[fa.NX, fa.NY + 1] + x[fa.NX + 1, fa.NY]);
			}
		} // set_bnd()
	}

	[DebuggerDisplay("{m_strName}")]
	public class Matrix
	{
		public string m_strName = "";
		public float[,] array = new float[fa.NX+2,fa.NY+2];

		public Matrix(string ip_strName)
		{
			m_strName = ip_strName;
		}

		public float this[int i, int j]
		{
			get
			{
				return array[i, j];
			}
			set
			{
				array[i,j] = value;
			}
		}

		public void Copy(Matrix ip_MatrixDest)
		{
			Array.Copy(array, ip_MatrixDest.array, ip_MatrixDest.array.Length);
		}
	}

	public class DrawPanel : Panel
	{
		public DrawPanel()
		{
			DoubleBuffered = true;
		}

		protected override void OnPaintBackground(PaintEventArgs pevent)
		{
		}
	}
}