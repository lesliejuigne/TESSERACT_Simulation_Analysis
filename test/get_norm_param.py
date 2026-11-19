import uproot
import pandas as pd
import warnings
from collections import defaultdict

warnings.simplefilter("ignore")

class GetNormParam:
    def __init__(self, root_path: str, objname: str):
        """
        Load a TMacro-like object from a ROOT file and extract the material masses and the initial 
        beamOn number.

        - objname determines the output:
            * 'runMacro' → extract beamOn number
            * 'geometryTable' → compute total mass per material
            
        How to use:
        >>> gnp = GetNormParam("path/to/file.root", "runMacro")
        >>> print(gnp.beamon_number)
        >>> gnp = GetNormParam("path/to/file.root", "geometryTable")
        >>> print(gnp.total_mass)
        """
        
        self.root_path = root_path
        self.objname = objname

        # Load table
        self.df = self._load_geometry_uproot()

        # Compute depending on object type
        self.beamon_number = None
        self.total_mass = None

        if self.objname == "runMacro":
            self.beamon_number = self._extract_beamon_number()
        elif self.objname == "geometryTable":
            self.total_mass = self._compute_total_mass()

    def _load_geometry_uproot(self):
        """load a TMacro into a pandas DataFrame."""
        f = uproot.open(self.root_path)

        if self.objname not in f:
            raise KeyError(f"'{self.objname}' not found in file. Run f.classnames() to locate it.")
        obj = f[self.objname]

        if "fLines" not in obj.all_members:
            raise KeyError("object has no member 'fLines' — not a TMacro-like object?")
        lines_list = obj.member("fLines")

        # Build list of strings
        lines = []
        for item in lines_list:
            if isinstance(item, (bytes, bytearray)):
                s = item.decode("utf-8", errors="ignore")
            else:
                s = str(item)
            if s is not None:
                lines.append(s.rstrip("\n"))

        rows = [ln.strip().split() for ln in lines if ln.strip()]
        df = pd.DataFrame(rows)
        
        return df

    def _extract_beamon_number(self):
        """Extract the number after '/run/beamOn' from the DataFrame."""
        mask = self.df.iloc[:,0] == '/run/beamOn'
        if mask.any():
            val = self.df.loc[mask, self.df.columns[1]].iloc[0]
            try:
                return int(val)
            except (ValueError, TypeError):
                return None
        return None

    def _compute_total_mass(self):
        """Compute the total mass per material from the DataFrame."""
        if self.df.empty:
            return {}

        material_col = self.df.columns[-1]
        mass_col = self.df.columns[1]

        total_mass = defaultdict(float)
        for _, row in self.df.iterrows():
            try:
                mass = float(row[mass_col])
                material = str(row[material_col])
                total_mass[material] += mass
            except (ValueError, TypeError):
                continue

        return dict(total_mass)
