/**********************************************************************
  tensorviswidget.h Base class for crystal builder dockwidgets

  Copyright (C) 2011 by David C. Lonie

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.openmolecules.net/>

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
 ***********************************************************************/

#include "tensorviswidget.h"

#include <avogadro/glwidget.h>
#include <avogadro/primitivelist.h>
#include <avogadro/atom.h>

#include <QtGui/QMessageBox>

#include <stdio.h>

namespace TensorVis
{

  TensorVisWidget::TensorVisWidget()
    : QDockWidget()
  {
    ui.setupUi(this);

    connect(ui.push_render, SIGNAL(clicked()), SIGNAL(generateMesh()));
    connect(ui.push_remove, SIGNAL(clicked()), SIGNAL(clearMesh()));
    connect(ui.push_selectAtom, SIGNAL(clicked()),
            SLOT(useSelectedAtomAsOrigin()));

    ui.edit_tensor->setCurrentFont(QFont("Monospace"));

    Eigen::Matrix3f tensor;
    tensor <<
      88.921500, 46.045900, -0.777100,
      17.744000, 18.222400, -0.248500,
      -0.900100, -0.432300, -81.653700;
    this->ui.edit_tensor->setMatrix(tensor.cast<double>());

    refreshGui();
  }

  TensorVisWidget::~TensorVisWidget()
  {
  }

  Eigen::Matrix3f TensorVisWidget::tensor()
  {
    return this->ui.edit_tensor->matrix().cast<float>();
  }

  Eigen::Vector3f TensorVisWidget::origin()
  {
    return Eigen::Vector3f(this->ui.spin_x->value(),
                           this->ui.spin_y->value(),
                           this->ui.spin_z->value());
  }

  Avogadro::Color3f TensorVisWidget::posColor()
  {
    return Avogadro::Color3f(0.2f, 0.2f, 0.7f);
  }

  Avogadro::Color3f TensorVisWidget::negColor()
  {
    return Avogadro::Color3f(0.7f, 0.2f, 0.2f);
  }

  void TensorVisWidget::refreshGui()
  {
    this->ui.edit_tensor->resetMatrix();
    this->ui.edit_tensor->setCurrentFont(QFont("Monospace"));
  }

  void TensorVisWidget::useSelectedAtomAsOrigin()
  {
    Avogadro::GLWidget *gl = Avogadro::GLWidget::current();
    Avogadro::PrimitiveList atoms =
        gl->selectedPrimitives().subList(Avogadro::Primitive::AtomType);

    if (atoms.count() != 1) {
      QMessageBox::information(gl, "Not so fast...",
                               "Please select one (and only one) atom to use "
                               "for the tensor origin.");
      return;
    }

    // Static cast is ok because of the above sublist call
    Avogadro::Atom *atom = static_cast<Avogadro::Atom*>(atoms.list().first());

    this->ui.spin_x->setValue(atom->pos()->x());
    this->ui.spin_y->setValue(atom->pos()->y());
    this->ui.spin_z->setValue(atom->pos()->z());
  }


}

#include "tensorviswidget.moc"
