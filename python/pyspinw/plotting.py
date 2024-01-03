import numpy as np
from vispy import scene
from vispy.color import color_array
from itertools import chain
from vispy.visuals.filters import ShadingFilter, WireframeFilter
from vispy.geometry import create_sphere
import copy
from scipy.spatial.transform import Rotation


class SuperCellSimple:
    def __init__(self, swobj, extent=(1,1,1), plot_mag=True, plot_bonds=False, plot_atoms=True, plot_cell=True, plot_axes=True, plot_plane=True, ion_type=None):
        # init with sw obj - could get NExt from object if not explicitly provide (i.e. make default None)
        self.do_plot_mag=plot_mag
        self.do_plot_bonds=plot_bonds
        self.do_plot_atoms=plot_atoms
        self.do_plot_cell=plot_cell
        self.do_plot_axes=plot_axes
        self.do_plot_plane=plot_plane
        self.do_plot_ion = ion_type is not None
        self.ion_type = ion_type  # "aniso" or "g"

        # get properties from swobj
        unit_cell = UnitCellSimple(basis_vec=swobj.basisvector().T)
        # add atoms
        _, single_ion = swobj.intmatrix('plotmode', True,'extend',False,'sortDM',False,'zeroC',False,'nExt',[1, 1, 1])
        for iatom, atom_idx in enumerate(swobj.atom()['idx'].flatten().astype(int)):
            imatom = np.flatnonzero(swobj.matom()['idx'].flatten().astype(int)== atom_idx)
            if len(imatom) > 0:
                S = swobj.matom()['S'].flat[imatom]
            else:
                S = None
            color = swobj.unit_cell['color'][:,atom_idx-1]/255
            label = swobj.unit_cell['label'][atom_idx-1]
            unit_cell.add_atom(AtomSimple(swobj.atom()['r'][:,iatom], S=S, size=0.35, color=color, label=label,
                                          gtensor_mat=single_ion['g'][:,:,iatom], aniso_mat=single_ion['aniso'][:,:,iatom]))
            
        # add bonds - only plot bonds for which there is a mat_idx
        bond_idx = np.squeeze(swobj.coupling['idx'])
        for ibond in np.unique(bond_idx[np.any(swobj.coupling['mat_idx'], axis=0)]):
            i_dl = np.squeeze(bond_idx==ibond)
            mat_idx = swobj.coupling['mat_idx'][0, np.argmax(i_dl)] - 1
            unit_cell.add_bond_vertices(ibond, np.squeeze(swobj.coupling['atom1'])[i_dl]-1, np.squeeze(swobj.coupling['atom2'])[i_dl]-1, swobj.coupling['dl'].T[i_dl], color=swobj.matrix['color'][:,mat_idx]/255)
        
        # generate unit cells in supercell
        self.unit_cells = []
        self.extent = np.asarray(extent)
        self.int_extent = np.ceil(self.extent).astype(int) + 1  # to plot additional unit cell along each dimension to get atoms on cell boundary
        for zcen in range(self.int_extent[2]):
            for ycen in range(self.int_extent[1]):
                for xcen in range(self.int_extent[0]):
                    self.unit_cells.append(unit_cell.translate_by([xcen, ycen, zcen]))
                    
        # get magnetic structure for all spins in supercell
        self.n = np.zeros(3)
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
        self.rotation_plane_radius = 0.3*self.cell_scale_abc_to_xyz
        self.ion_radius = 0.3*self.cell_scale_abc_to_xyz

    def transform_points_abc_to_xyz(self, points):
        return points @ self.basis_vec
        
    def apply_magnetic_structure(self, swobj):
        magstr = swobj.magstr('NExt', [int(ext) for ext in self.int_extent])
        mj = magstr['S']
        if not np.any(mj):
            print('No magnetic structure defined')
            self.do_plot_mag = False
            self.do_plot_plane = False
            return
        icol = 0
        for uc in self.unit_cells:
            # think is generated in correct order...
            for atom in uc.atoms:
                if atom.S is not None:
                    atom.moment = mj[:,icol]
                    icol+=1
        self.n = np.asarray(magstr['n'])  # plane of rotation of moment

    def plot(self):
        canvas = scene.SceneCanvas(bgcolor='white',  show=True) # size=(600, 600),
        view = canvas.central_widget.add_view()
        view.camera = scene.cameras.TurntableCamera() # fov=5)
        
        pos, colors, sizes, labels, iremove = self.get_atomic_positions_xyz_in_supercell()
        subvisuals = []
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
        if self.do_plot_plane:
            subvisuals.extend(self.make_rotation_plane_visuals())
        if self.do_plot_ion:
            subvisuals.extend(self.make_ion_visuals())
        # display
        if subvisuals:
            scene.Compound(subvisuals=subvisuals, parent=view.scene)
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
                    arrows=arrows, arrow_type='angle_60', arrow_size=3.,
                    width=3., antialias=False, arrow_color=arrow_color,
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
        pos = np.array([atom.pos + cell.origin for cell in self.unit_cells for atom in cell.atoms])
        sizes = np.array([atom.size for cell in self.unit_cells for atom in cell.atoms])
        colors = np.array([atom.color for cell in self.unit_cells for atom in cell.atoms]).reshape(-1,3)
        # remove points and vertices beyond extent
        pos, iremove = self._remove_points_outside_extent(pos)
        sizes = np.delete(sizes, iremove)
        colors = np.delete(colors, iremove, axis=0)
        # transfrom to xyz
        pos = self.transform_points_abc_to_xyz(pos)
        # get atomic labels
        labels = [atom.label for cell in self.unit_cells for atom in cell.atoms]
        labels = np.delete(labels, iremove).tolist()
        return pos, colors, sizes, labels, iremove
    
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
    
    def make_rotation_plane_visuals(self, npts=15):
        # generate vertices of disc with normal // [0,0,1]
        theta = np.linspace(0, 2*np.pi,npts-1)
        disc_verts = np.zeros((npts, 3))
        disc_faces = np.zeros((npts-2, 3), dtype=int)
        disc_verts[1:,0] = self.rotation_plane_radius*np.cos(theta)
        disc_verts[1:,1] = self.rotation_plane_radius*np.sin(theta)
        # rotate given normal
        rot_mat = get_rotation_matrix(self.n)
        disc_verts = rot_mat.dot(disc_verts.T).T
        # label faces
        disc_faces[:,1] = np.arange(1,npts-1)
        disc_faces[:,2] = np.arange(2,npts)
        disc_visuals = []
        for cell in self.unit_cells:
            for atom in cell.atoms:
                if atom.S is not None and atom.S > 0:
                    centre = (atom.pos + cell.origin).reshape(1,-1) # i.e. make 2D array
                    centre, _ = self._remove_points_outside_extent(centre)
                    if centre.size > 0:
                        # atom in extents
                        centre = self.transform_points_abc_to_xyz(centre)
                        mesh = scene.visuals.Mesh(vertices=disc_verts + centre, faces=disc_faces, color=color_array.Color(color=atom.color, alpha=0.25))
                        disc_visuals.append(mesh)
        return disc_visuals
    
    def make_ion_visuals(self, npts=7):
        # get mesh for a sphere
        meshdata = create_sphere(radius=self.ion_radius, rows=npts, cols=npts)
        verts = meshdata.get_vertices()
        faces = meshdata.get_faces()
        ellips_visuals = []
        for cell in self.unit_cells:
            for atom in cell.atoms:
                    centre = (atom.pos + cell.origin).reshape(1,-1) # i.e. make 2D array
                    centre, _ = self._remove_points_outside_extent(centre)
                    if centre.size > 0:
                        # atom in extents
                        transform = atom.get_transform(tensor=self.ion_type)
                        if np.any(transform>0):
                            centre = self.transform_points_abc_to_xyz(centre)
                            this_verts = verts @ transform + centre
                            mesh = scene.visuals.Mesh(vertices=this_verts, faces=faces, color=color_array.Color(color=atom.color, alpha=0.25))
                            wireframe_filter = WireframeFilter(color=3*[0.7])
                            mesh.attach(wireframe_filter)
                            ellips_visuals.append(mesh)
        return ellips_visuals

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
        self.origin = np.array(origin)
        self.basis_vec = basis_vec  # each col is a basis vector
    
    def add_atom(self, atom):
        self.atoms.append(atom)

    def add_bond_vertices(self, name, atom1_idx, atom2_idx, dl, color):
        # store tuple of vertices in same way as spins returned
        self.bonds[name] = {'verts': np.array([(self.atoms[atom1_idx[ibond]].pos, 
                                                self.atoms[atom2_idx[ibond]].pos + dl) for ibond, dl in enumerate(np.asarray(dl))]).reshape(-1,3)}
        self.bonds[name]['color'] = color

    def translate_by(self, origin):
        return UnitCellSimple(copy.deepcopy(self.atoms), self.bonds, origin, basis_vec=self.basis_vec)

    def get_bond_vertices(self, bond_name):
        return self.bonds[bond_name]['verts'] + self.origin

    def get_bond_color(self, bond_name):
        return self.bonds[bond_name]['color']


class AtomSimple:
    def __init__(self, position, S=None, moment=np.zeros(3), size=0.2, gtensor_mat=None, aniso_mat=None, label='atom', color="blue"):
        self.pos = np.asarray(position)
        self.S = S
        self.moment=moment
        self.gtensor = gtensor_mat
        self.aniso = aniso_mat
        self.size = size
        self.color = color
        self.label = label
        self.spin_scale = 0.3
        
    def get_transform(self, tensor='aniso'):
        if tensor=="aniso":
            mat = self.aniso
        else:
            mat = self.gtensor
        # diagonalise so can normalise eigenvalues 
        evals, evecs = np.linalg.eig(mat)
        if tensor=="aniso":
            # take inverse of eigenvals as large number should produce a samll axis
            evals = 1/evals
        # scale such that max eigval is 1
        evals = evals/np.max(abs(evals))
        return evecs @ np.diag(evals) @ np.linalg.inv(evecs)


def get_rotation_matrix(vec2, vec1=np.array([0,0,1])):
    vec1 = vec1/np.linalg.norm(vec1)  # unit vectors
    vec2 = vec2/np.linalg.norm(vec2)
    if np.arccos(np.clip(np.dot(vec1.flat, vec2.flat), -1.0, 1.0)) > 1e-5:
        r = Rotation.align_vectors(vec2.reshape(1,-1), vec1.reshape(1,-1))  # matrix to rotate vec1 onto vec2
        return r[0].as_matrix()
    else:
        # too small a difference for above algorithm, just return identity
        return np.eye(3)