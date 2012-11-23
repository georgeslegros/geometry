using System;
using System.Windows;

namespace Foo
{
    public class Geometry
    {
        //calculates the theta value used to calculate the endpoints on the circle to render lines and arcs.
        public static double CalculateTheta(double sectionCount)
        {
            return 2 * Math.PI / sectionCount;
        }

        public static double Phi(double theta)
        {
            return -(Math.PI / 4.0) - (theta / 2.0);
        }

        public static Point CalculateCircleEndPoint(double angle, double radius, Point circleCenter)
        {
            double x = (radius * (Math.Cos(angle))) + circleCenter.X;
            double y = (radius * (Math.Sin(angle))) + circleCenter.Y;
            return new Point(x, y);
        }

        // Find the points where the two circles intersect.
        public static int FindCircleCircleIntersections(Point p1, double radius0, Point p2, double radius1,
                                                        out Point intersection1, out Point intersection2)
        {
            // Find the distance between the centers.
            var dx = p1.X - p2.X;
            var dy = p1.Y - p2.Y;
            double dist = Math.Sqrt(dx * dx + dy * dy);

            // See how manhym solutions there are.
            if (dist > radius0 + radius1)
            {
                // No solutions, the circles are too far apart.
                intersection1 = new Point(float.NaN, float.NaN);
                intersection2 = new Point(float.NaN, float.NaN);
                return 0;
            }
            if (dist < Math.Abs(radius0 - radius1))
            {
                // No solutions, one circle contains the other.
                intersection1 = new Point(float.NaN, float.NaN);
                intersection2 = new Point(float.NaN, float.NaN);
                return 0;
            }
            if ((dist == 0) && (radius0 == radius1))
            {
                // No solutions, the circles coincide.
                intersection1 = new Point(float.NaN, float.NaN);
                intersection2 = new Point(float.NaN, float.NaN);
                return 0;
            }

            // Find a and h.
            double a = (radius0 * radius0 -
                        radius1 * radius1 + dist * dist) / (2 * dist);
            double h = Math.Sqrt(radius0 * radius0 - a * a);

            // Find P2.
            double cx2 = p1.X + a * (p2.X - p1.X) / dist;
            double cy2 = p1.Y + a * (p2.Y - p1.Y) / dist;

            // Get the points P3.
            intersection1 = new Point(
                (float)(cx2 + h * (p2.Y - p1.Y) / dist),
                (float)(cy2 - h * (p2.X - p1.X) / dist));
            intersection2 = new Point(
                (float)(cx2 - h * (p2.Y - p1.Y) / dist),
                (float)(cy2 + h * (p2.X - p1.X) / dist));

            // See if we have 1 or 2 solutions.
            if (dist == radius0 + radius1) return 1;
            return 2;
        }

        // Find the points of intersection.
        public static int FindLineCircleIntersections(Point circleCenter, double radius, Point point1, Point point2,
                                                      out Point intersection1, out Point intersection2)
        {
            double dx, dy, A, B, C, det, t;

            dx = point2.X - point1.X;
            dy = point2.Y - point1.Y;

            A = dx * dx + dy * dy;
            B = 2 * (dx * (point1.X - circleCenter.X) + dy * (point1.Y - circleCenter.Y));
            C = (point1.X - circleCenter.X) * (point1.X - circleCenter.X) +
                (point1.Y - circleCenter.Y) * (point1.Y - circleCenter.Y) - radius * radius;

            det = B * B - 4 * A * C;
            if ((A <= 0.0000001) || (det < 0))
            {
                // No real solutions.
                intersection1 = new Point(float.NaN, float.NaN);
                intersection2 = new Point(float.NaN, float.NaN);
                return 0;
            }
            else if (det == 0)
            {
                // One solution.
                t = -B / (2 * A);
                intersection1 = new Point(point1.X + t * dx, point1.Y + t * dy);
                intersection2 = new Point(float.NaN, float.NaN);
                return 1;
            }
            else
            {
                // Two solutions.
                t = (float)((-B + Math.Sqrt(det)) / (2 * A));
                intersection1 = new Point(point1.X + t * dx, point1.Y + t * dy);
                t = (float)((-B - Math.Sqrt(det)) / (2 * A));
                intersection2 = new Point(point1.X + t * dx, point1.Y + t * dy);
                return 2;
            }
        }

        public static double GetDistanceBetweenTwoPoints(Point point1, Point point2)
        {
            double dx = point1.X - point2.X;
            double dy = point1.Y - point2.Y;
            double dist = Math.Sqrt(dx * dx + dy * dy);
            return dist;
        }

        public static Point GetCircleCenterGivenTwoPointsAndRadius(Point a, Point b, double radius)
        {
            double d = GetDistanceBetweenTwoPoints(a, b);

            double factor = Math.Pow(radius, 2) - (Math.Pow(d, 2) / 4);

            double x = ((b.X - a.X) / 2) + Math.Sqrt(factor / (1 + Math.Pow((b.X - a.X) / (b.Y - a.Y), 2)));
            double y = ((b.Y - a.Y) / 2) + Math.Sqrt(factor / (1 + Math.Pow((b.Y - a.Y) / (b.X - a.X), 2)));

            return new Point(x, y);
        }

        public static Point MovePointAlongLine(Point point, double angle, double distance)
        {
            var radAngle = angle * Math.PI / 180.0;
            var x = point.X + Math.Cos(radAngle) * distance;
            var y = point.Y + Math.Sin(radAngle) * distance;
            return new Point(x, y);
        }

        public static Point MovePointAlongCircle(Point point, Point centerOfCircle, double angleToMove)
        {
            //  Get the point in polar coords.
            double radius, phi;
            CartesianToPolar(point, centerOfCircle, out radius, out phi);

            //  Adjust the angle.
            phi += angleToMove;

            //  Back to Cartesian
            double x = (float)Math.Cos(phi) * radius;
            double y = (float)Math.Sin(phi) * radius;

            //  Offset.
            return new Point(x, y).Add(centerOfCircle);
        }

        public static Point GetArcMiddle(Point start, Point end, Point center)
        {
            double startPolarRadius, startPolarAngle;
            CartesianToPolar(start, center, out startPolarRadius, out startPolarAngle);

            double endPolarRadius, endPolarAngle;
            CartesianToPolar(end, center, out endPolarRadius, out endPolarAngle);

            if (startPolarAngle > endPolarAngle)
                endPolarAngle += 2.0 * Math.PI;

            double radius = (startPolarRadius + endPolarRadius) / 2.0;
            double phi = (startPolarAngle + endPolarAngle) / 2.0;
            
            return PolarToCartesian(center, radius, phi);
        }

        public static void CartesianToPolar(Point point, Point center, out double radius, out double phi)
        {
            var p = point.Minus(center);
            radius = Math.Sqrt((p.X * p.X) + (p.Y * p.Y));
            phi = Math.Atan2(p.Y, p.X);
        }

        public static Point PolarToCartesian(Point center, double radius, double phi)
        {
            var x = radius * Math.Cos(phi);
            var y = radius * Math.Sin(phi);
            return new Point(x, y).Add(center);
        }

        public static double NormalizeAngle(double angle)
        {
            double newAngle = angle;
            while (newAngle <= 0) newAngle += 360;
            while (newAngle > 360) newAngle -= 360;
            return newAngle;
        }
    }
}