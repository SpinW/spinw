import numpy as np
from vispy import scene
from vispy.color import color_array
from itertools import chain
from vispy.visuals.filters import ShadingFilter, WireframeFilter
from vispy.geometry import create_sphere
import copy


class SuperCellSimple:
    def __init__(self, swobj, extent=(1,1,1), plot_mag=True, plot_bonds=False, plot_atoms=True, plot_cell=True, plot_axes=True):
        # init with sw obj - could get NExt from object if not explicitly provide (i.e. make default None)
        self.do_plot_mag=plot_mag
        self.do_plot_bonds=plot_bonds
        self.do_plot_atoms=plot_atoms
        self.do_plot_cell=plot_cell
        self.do_plot_axes=plot_axes
        
        # get properties from swobj
        unit_cell = UnitCellSimple(basis_vec=swobj.basisvector().T)
        # add atoms
        for iatom, atom_idx in enumerate(swobj.atom()['idx'].flatten().astype(int)):
            imatom = np.flatnonzero(swobj.matom()['idx'].flatten().astype(int)== atom_idx)
            if len(imatom) > 0:
                S = swobj.matom()['S'].flat[imatom]
            else:
                S = None
            # get color
            unit_cell.add_atom(AtomSimple(swobj.atom()['r'][:,iatom], S=S, size=0.35, color=swobj.unit_cell['color'][:,iatom], label=swobj.unit_cell['label'][iatom]))
        # only plot bonds for which there is a mat_idx
        bond_idx = np.squeeze(swobj.coupling['idx'])
        for ibond in np.unique(bond_idx[np.any(swobj.coupling['mat_idx'], axis=0)]):
            i_dl = np.squeeze(bond_idx==ibond)
            mat_idx = swobj.coupling['mat_idx'][0, np.argmax(i_dl)] - 1
            unit_cell.add_bond_vertices(ibond, np.squeeze(swobj.coupling['atom1'])[i_dl]-1, np.squeeze(swobj.coupling['atom2'])[i_dl]-1, swobj.coupling['dl'].T[i_dl], color=swobj.matrix['color'][:,mat_idx]/365)
        
        # generate unit cells in supercell
        self.unit_cells = []
        self.extent = np.asarray(extent)
        self.int_extent = np.ceil(self.extent).astype(int) + 1  # to plot additional unit cell along each dimension to get atoms on cell boundary
        for zcen in range(self.int_extent[2]):
            for ycen in range(self.int_extent[1]):
                for xcen in range(self.int_extent[0]):
                    self.unit_cells.append(unit_cell.translate_by([xcen, ycen, zcen]))
                    
        # get magnetic structure for all spins in supercell
        self.apply_magnetic_structure(swobj)
        
        # transforms
        self.basis_vec = unit_cell.basis_vec
                
        # scale factors
        self.abc = np.sqrt(np.sum(self.basis_vec**2, axis=1))
        self.cell_scale_abc_to_xyz = min(self.abc)
        self.supercell_scale_abc_to_xyz = min(self.abc*self.extent)
        self.bond_width = 5
        self.spin_scale = 1
        self.arrow_width = 8
        self.arrow_head_size = 5
        self.font_size = 30
        self.axes_font_size = 50
        self.atom_alpha = 0.5

    def transform_points_abc_to_xyz(self, points):
        return points @ self.basis_vec
        
    def apply_magnetic_structure(self, swobj):
        mj = swobj.magstr('NExt', [int(ext) for ext in self.int_extent])['S']
        if not np.any(mj):
            print('No magnetic structure defined')
            self.plot_mag = False
            return
        icol = 0
        for uc in self.unit_cells:
            # think is generated in correct order...
            for atom in uc.atoms:
                if atom.S is not None:
                    atom.moment = mj[:,icol]
                    icol+=1

    def plot(self):
        canvas = scene.SceneCanvas(bgcolor='white',  show=True) # size=(600, 600),
        view = canvas.central_widget.add_view()
        view.camera = scene.cameras.TurntableCamera() # fov=5)
        
        pos, colors, sizes, labels, iremove = self.get_atomic_positions_xyz_in_supercell()
        if self.do_plot_cell:
            self.plot_unit_cell_box(view.scene)  # plot girdlines for unit cell boundaries
        if self.do_plot_mag:
            self.plot_magnetic_structure(view.scene, pos, colors, iremove)
        if self.do_plot_atoms:
            self.plot_atoms(view.scene, pos, colors, sizes, labels)
        if self.do_plot_bonds:
            self.plot_bonds(view.scene)
        if self.do_plot_axes:
            self.plot_cartesian_axes(view.scene)
        # display
        view.camera.set_range()  # centers camera on middle of data and auto-scales extent
        canvas.app.run()
        return canvas, view.scene

    def plot_cartesian_axes(self, canvas_scene):
        pos = np.array([[0., 0., 0.], [1., 0., 0.],
                        [0., 0., 0.], [0., 1., 0.],
                        [0., 0., 0.], [0., 0., 1.],
                        ])*0.5
        pos = pos - 0.5*np.ones(3)
        pos = self.transform_points_abc_to_xyz(pos)
        arrows = np.c_[pos[0::2], pos[1::2]]

        line_color = ['red', 'red', 'green', 'green', 'blue', 'blue']
        arrow_color = ['red', 'green', 'blue']

        scene.visuals.Arrow(pos=pos, parent=canvas_scene, connect='segments',
                    arrows=arrows, arrow_type='triangle_60', arrow_size=3.,
                    width=3., antialias=True, arrow_color=arrow_color,
                    color=line_color)
        scene.visuals.Text(pos=self.transform_points_abc_to_xyz(0.7*np.eye(3)-0.5), parent=canvas_scene, text=["a", "b", "c"], color=arrow_color, 
                           font_size=self.axes_font_size*self.supercell_scale_abc_to_xyz)


    def plot_unit_cell_box(self, canvas_scene):
        for zcen in range(self.int_extent[2]):
            for ycen in range(self.int_extent[1]):
                    scene.visuals.Line(pos = self.transform_points_abc_to_xyz(np.array([[0, ycen, zcen], [np.ceil(self.extent[0]), ycen, zcen]])),
                                       parent=canvas_scene, color=color_array.Color(color="k", alpha=0.25)) # , method="gl")
        for xcen in range(self.int_extent[0]):
            for ycen in range(self.int_extent[1]):
                    scene.visuals.Line(pos = self.transform_points_abc_to_xyz(np.array([[xcen, ycen, 0], [xcen, ycen, np.ceil(self.extent[2])]])),
                                       parent=canvas_scene, color=color_array.Color(color="k", alpha=0.25)) # , method="gl")
        for xcen in range(self.int_extent[0]):
            for zcen in range(self.int_extent[2]):
                    scene.visuals.Line(pos = self.transform_points_abc_to_xyz(np.array([[xcen, 0, zcen], [xcen, np.ceil(self.extent[1]), zcen]])),
                                       parent=canvas_scene, color=color_array.Color(color="k", alpha=0.25)) # , method="gl")

    def get_atomic_positions_xyz_in_supercell(self):
        pos = np.array([atom._pos + cell._origin for cell in self.unit_cells for atom in cell.atoms])
        sizes = np.array([atom._size for cell in self.unit_cells for atom in cell.atoms])
        colors = np.array([atom._color for cell in self.unit_cells for atom in cell.atoms]).reshape(-1,3)
        # remove points and vertices beyond extent
        pos, iremove = self._remove_points_outside_extent(pos)
        sizes = np.delete(sizes, iremove)
        colors = np.delete(colors, iremove, axis=0)
        # transfrom to xyz
        pos = self.transform_points_abc_to_xyz(pos)
        # get atomic labels
        labels = [atom.label for cell in self.unit_cells for atom in cell.atoms]
        labels = np.delete(labels, iremove).tolist()
        return pos, colors/365, sizes, labels, iremove
        
    def plot_magnetic_structure(self, canvas_scene, pos, colors, iremove):
        mj = np.array([atom.moment for cell in self.unit_cells for atom in cell.atoms]) # already in xyz
        mj = self.spin_scale*np.delete(mj, iremove, axis=0)
        verts = np.c_[pos, pos+mj]  # natom x 6
        scene.visuals.Arrow(pos=verts.reshape(-1,3), parent=canvas_scene, connect='segments',
            arrows=verts, arrow_size=self.arrow_head_size,
            width=self.arrow_width, antialias=True, 
            arrow_type='triangle_60',
            color=np.repeat(colors, 2, axis=0).tolist(),
            arrow_color= colors.tolist())

    def plot_atoms(self, canvas_scene, pos, colors, sizes, labels):
        scene.visuals.Markers(
                    pos=pos,
                    size=sizes*self.cell_scale_abc_to_xyz,
                    antialias=0,
                    face_color= colors,
                    edge_color='white',
                    edge_width=0,
                    scaling=True,
                    spherical=True,
                    alpha=self.atom_alpha,
                    parent=canvas_scene)
        # labels
        scene.visuals.Text(pos=pos, parent=canvas_scene, text=labels, color="white", font_size=self.font_size*self.cell_scale_abc_to_xyz)

    def plot_bonds(self, canvas_scene):
        for bond_name in self.unit_cells[0].bonds:
            color = self.unit_cells[0].get_bond_color(bond_name)
            verts = np.array([unit_cell.get_bond_vertices(bond_name) for unit_cell in self.unit_cells]).reshape(-1,3)
            verts, _ = self._remove_vertices_outside_extent(verts)
            verts = self.transform_points_abc_to_xyz(verts)
            scene.visuals.Line(pos=verts, parent=canvas_scene, connect='segments', 
                               width=self.bond_width, color=color)

    def _remove_vertices_outside_extent(self, verts):
        # DO THIS BEFORE CONVERTING TO XYZ
        # remove pairs of verts that correpsond to a line outside the extent
        _, iremove = self._remove_points_outside_extent(verts)
        iatom2 = (iremove % 2).astype(bool) # end point of pair of vertices
        # for atom2 type vertex we need to remove previous row (atom1 vertex)
        # for atom1 type vertex we need to remove the subsequent row (atom2 vertex)
        iremove = np.hstack((iremove, iremove[iatom2]-1, iremove[~iatom2]+1))
        return np.delete(verts, iremove, axis=0), iremove

    def _remove_points_outside_extent(self, points):
        # DO THIS BEFORE CONVERTING TO XYZ
        iremove = np.flatnonzero(np.logical_or(np.any(points<0, axis=1), np.any((points - self.extent)>0, axis=1)))
        return np.delete(points, iremove, axis=0), iremove

class UnitCellSimple:
    def __init__(self, atoms_list=[], bonds={}, origin=np.array([0,0,0]), basis_vec=np.eye(3)):
        self.atoms = atoms_list
        self.bonds = bonds
        self._origin = np.array(origin)
        self.basis_vec = basis_vec  # each col is a basis vector
    
    def add_atom(self, atom):
        self.atoms.append(atom)

    def add_bond_vertices(self, name, atom1_idx, atom2_idx, dl, color):
        # store tuple of vertices in same way as spins returned
        self.bonds[name] = {'verts': np.array([(self.atoms[atom1_idx[ibond]]._pos, 
                                                self.atoms[atom2_idx[ibond]]._pos + dl) for ibond, dl in enumerate(np.asarray(dl))]).reshape(-1,3)}
        self.bonds[name]['color'] = color

    def translate_by(self, origin):
        return UnitCellSimple(copy.deepcopy(self.atoms), self.bonds, origin, basis_vec=self.basis_vec)

    def get_bond_vertices(self, bond_name):
        return self.bonds[bond_name]['verts'] + self._origin

    def get_bond_color(self, bond_name):
        return self.bonds[bond_name]['color']


class AtomSimple:
    def __init__(self, position, S=None, moment=np.zeros(3), size=0.2, gtensor_mat=None, aniso_mat=None, n=None, label='atom', color="blue"):
        self._pos = np.asarray(position)
        self.S = S
        self.moment=moment
        self._gtensor = gtensor_mat
        self._aniso = aniso_mat
        self._n = np.asarray(n)
        self._size = size
        self._color = color
        self.label = label
        self.spin_scale = 0.3
        
    def get_transform(self, tensor='aniso'):
        if tensor=="aniso":
            return self._aniso
        else:
            return self._gtensor
