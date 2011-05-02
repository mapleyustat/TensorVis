/**********************************************************************
  TensorMeshEngine - Engine for display meshes build from tensors

  Copyright (C) 2011 by David C. Lonie

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 **********************************************************************/

#include "tensormeshengine.h"

#include <avogadro/color.h>
#include <avogadro/mesh.h>
#include <avogadro/molecule.h>
#include <avogadro/painter.h>
#include <avogadro/painterdevice.h>

#include <QDebug>

using namespace Avogadro;

namespace TensorVis {

  TensorMeshEngine::TensorMeshEngine(QObject *parent) : Engine(parent)
  {
  }

  Engine *TensorMeshEngine::clone() const
  {
    TensorMeshEngine *engine = new TensorMeshEngine(parent());
    engine->setAlias(alias());
    engine->setEnabled(isEnabled());

    return engine;
  }

  TensorMeshEngine::~TensorMeshEngine()
  {
  }

  bool TensorMeshEngine::renderOpaque(PainterDevice *pd)
  {
    Painter *p = pd->painter();
    const Molecule *mol = pd->molecule();
    QList<Mesh*> meshes = mol->meshes();

    // Find and render all tensor meshes
    for (QList<Mesh*>::const_iterator
           it = meshes.constBegin(),
           it_end = meshes.constEnd();
         it != it_end; ++it) {
      /// @todo this may break when translated
      if ((*it)->name().startsWith(tr("Tensor"))) {
        p->drawColorMesh(**it, 0);
      }
    }

    return true;
  }

  Engine::PrimitiveTypes TensorMeshEngine::primitiveTypes() const
  {
    return Engine::Surfaces;
  }

}

#include "tensormeshengine.moc"

Q_EXPORT_PLUGIN2(tensormeshengine, TensorVis::TensorMeshEngineFactory)
