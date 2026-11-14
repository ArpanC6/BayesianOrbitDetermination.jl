// Independent reference propagation with Orekit 13.1.9.
// Generates CSV states for the BayesianOrbitDetermination.jl propagator test.
//
// Usage: java -cp "lib/*:." OrekitReference <orekit-data-dir> <case: j2|j2drag> <out.csv>
//
// Case j2drag uses: Newtonian + J2 + exponential drag (spherical Earth,
// no flattening, rho0 = 2.789e-10 kg/m^3 at 200 km, H = 50 km, Cd = 2.2,
// area/mass = 0.01 m^2/kg), matching OrbitPropagatorOptions in the Julia package.

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Locale;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.gravity.J2OnlyPerturbation;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.models.earth.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

public class OrekitReference {

    public static void main(String[] args) throws Exception {
        if (args.length != 3) {
            System.err.println("usage: OrekitReference <orekit-data-dir> <j2|j2drag> <out.csv>");
            System.exit(2);
        }
        File dataDir = new File(args[0]);
        String which = args[1];
        File out = new File(args[2]);

        DataContext.getDefault().getDataProvidersManager()
                .addProvider(new DirectoryCrawler(dataDir));

        final double mu = 398600.4418e9;     // m^3/s^2
        final double rEq = 6378137.0;        // m
        final double j2 = 1.0826269e-3;

        Frame eme2000 = FramesFactory.getEME2000();
        AbsoluteDate epoch = new AbsoluteDate(2000, 1, 1, 12, 0, 0.0,
                TimeScalesFactory.getTAI());

        // KeplerianState(6900 km, 0.002, 52 deg, 30 deg, 40 deg, 0) -> ECI, km and km/s
        Vector3D r0 = new Vector3D(3205.830222953533, 4997.609661816432, 3488.022461182947).scalarMultiply(1000.0);
        Vector3D v0 = new Vector3D(-6.035351552351067, 0.6629103392316869, 4.597259758177123).scalarMultiply(1000.0);

        CartesianOrbit orbit = new CartesianOrbit(new PVCoordinates(r0, v0), eme2000, epoch, mu);

        DormandPrince853Integrator integrator =
                new DormandPrince853Integrator(1.0e-3, 600.0, 1.0e-6, 1.0e-13);
        NumericalPropagator propagator = new NumericalPropagator(integrator);
        propagator.setInitialState(new SpacecraftState(orbit, 1000.0)); // 1000 kg -> A = 10 m^2

        propagator.addForceModel(new NewtonianAttraction(mu));
        propagator.addForceModel(new J2OnlyPerturbation(mu, rEq, j2, eme2000));

        if (which.equals("j2drag")) {
            OneAxisEllipsoid earth = new OneAxisEllipsoid(rEq, 0.0,
                    FramesFactory.getITRF(IERSConventions.IERS_2010, false));
            SimpleExponentialAtmosphere atmosphere = new SimpleExponentialAtmosphere(
                    earth, 2.789e-10, 200000.0, 50000.0);
            propagator.addForceModel(new DragForce(atmosphere, new IsotropicDrag(10.0, 2.2)));
        } else if (!which.equals("j2")) {
            System.err.println("unknown case: " + which);
            System.exit(2);
        }

        try (PrintWriter w = new PrintWriter(new FileWriter(out))) {
            w.println("# Orekit 13.1.9"
                    + " case=" + which
                    + " integrator=DormandPrince853 absTol=1e-6 relTol=1e-13"
                    + " epoch=2000-01-01T12:00:00 TAI frame=EME2000");
            w.println("t_sec,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms");
            for (int t = 0; t <= 86400; t += 600) {
                SpacecraftState s = (t == 0)
                        ? propagator.getInitialState()
                        : propagator.propagate(epoch.shiftedBy((double) t));
                PVCoordinates pv = s.getPVCoordinates(eme2000);
                Vector3D p = pv.getPosition().scalarMultiply(1.0e-3);
                Vector3D v = pv.getVelocity().scalarMultiply(1.0e-3);
                w.println(String.format(Locale.ROOT,
                        "%d,%.12f,%.12f,%.12f,%.15f,%.15f,%.15f",
                        t, p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()));
            }
        }
        System.out.println("wrote " + out);
    }
}
