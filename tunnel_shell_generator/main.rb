  module Geometry
    def self.interpolate_points(p1, p2, t)
      Geom::Point3d.linear_combination(1.0 - t, p1, t, p2)
    end

    def self.project_point_to_segment(p, a, b)
      ab = b - a
      ab_len = ab.length
      return a if ab_len < 0.00001
      ab_unit = ab.normalize
      ap = p - a
      projection = ap.dot(ab_unit)
      if projection <= 0
        a
      elsif projection >= ab_len
        b
      else
        # ✅ اصلاح: ضرب صحیح Vector
        a + (ab_unit * projection)
      end
    end
  end
