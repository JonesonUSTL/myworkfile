#include "inp_parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fem {
namespace {

// 字符串去空白。
std::string trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

// 转大写，用于关键字无关大小写匹配。
std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return s;
}

// 逗号分隔解析。
std::vector<std::string> splitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(trim(item));
  return out;
}

// 关键字行结构。
struct KeywordLine {
  std::string key;
  std::unordered_map<std::string, std::string> params;
};

// 解析 *Keyword, a=b, c=d 形式。
KeywordLine parseKeyword(const std::string& line) {
  KeywordLine kw;
  auto parts = splitCSV(line);
  if (parts.empty()) return kw;
  kw.key = upper(parts[0]);
  for (size_t i = 1; i < parts.size(); ++i) {
    const auto p = parts[i];
    const auto eq = p.find('=');
    if (eq == std::string::npos) {
      kw.params[upper(p)] = "";
    } else {
      kw.params[upper(trim(p.substr(0, eq)))] = trim(p.substr(eq + 1));
    }
  }
  return kw;
}

// 单元类型映射。
ElementType parseElementType(const std::string& t) {
  const auto u = upper(t);
  if (u == "T3D2" || u == "T2D2") return ElementType::Truss2;
  if (u == "B31" || u == "B33") return ElementType::Beam2;
  if (u == "S4" || u == "S4R") return ElementType::Shell4;
  if (u == "C3D8" || u == "C3D8R") return ElementType::Solid8;
  throw std::runtime_error("Unsupported element type: " + t);
}

// 判断是否是纯数字（用于区分节点号与节点集名）。
bool isIntegerToken(const std::string& s) {
  if (s.empty()) return false;
  size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
  if (i >= s.size()) return false;
  for (; i < s.size(); ++i) {
    if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
  }
  return true;
}

// 在节点列表中查索引。
int nodeIndex(const Model& m, int id) {
  for (size_t i = 0; i < m.nodes.size(); ++i)
    if (m.nodes[i].id == id) return static_cast<int>(i);
  return -1;
}

// 向集合追加数据。
void appendSetItems(std::unordered_map<std::string, std::vector<int>>& sets, const std::string& name,
                    const std::vector<std::string>& p) {
  auto& vec = sets[name];
  for (const auto& it : p)
    if (!it.empty()) vec.push_back(std::stoi(it));
}

}  // namespace

Model parseInpFile(const std::string& filePath) {
  std::ifstream in(filePath);
  if (!in) throw std::runtime_error("Cannot open inp file: " + filePath);

  Model model;
  std::string line;
  std::string currentSection;
  std::string currentMaterial;
  std::string currentNset;
  std::string currentElset;
  std::string currentAmplitude;
  std::string currentPartName;
  std::string currentLoadAmplitude;
  std::string plasticMode = "BILINEAR";
  ElementType currentElementType = ElementType::Truss2;

  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line.rfind("**", 0) == 0) continue;

    if (line.front() == '*') {
      const auto kw = parseKeyword(line);
      currentSection = kw.key;

      if (kw.key == "*HEADING") {
        // 下一行若为普通文本，作为模型名称。
      } else if (kw.key == "*PART") {
        currentPartName = kw.params.count("NAME") ? kw.params.at("NAME") : "PART-UNNAMED";
        model.partNames.push_back(currentPartName);
      } else if (kw.key == "*END PART") {
        currentPartName.clear();
      } else if (kw.key == "*OUTPUT") {
        if (kw.params.count("FIELD") && kw.params.count("FREQUENCY")) {
          model.step.fieldOutputFrequency = std::max(1, std::stoi(kw.params.at("FREQUENCY")));
        }
      } else if (kw.key == "*FIELD OUTPUT") {
        if (kw.params.count("FREQUENCY")) model.step.fieldOutputFrequency = std::max(1, std::stoi(kw.params.at("FREQUENCY")));
      } else if (kw.key == "*MATERIAL") {
        currentMaterial = kw.params.count("NAME") ? kw.params.at("NAME") : "MAT1";
        model.materials[currentMaterial].name = currentMaterial;
      } else if (kw.key == "*ELEMENT") {
        currentElementType = kw.params.count("TYPE") ? parseElementType(kw.params.at("TYPE")) : ElementType::Truss2;
        currentElset = kw.params.count("ELSET") ? kw.params.at("ELSET") : "EALL";
      } else if (kw.key == "*NSET") {
        currentNset = kw.params.count("NSET") ? kw.params.at("NSET") : "NALL";
      } else if (kw.key == "*ELSET") {
        currentElset = kw.params.count("ELSET") ? kw.params.at("ELSET") : "EALL";
      } else if (kw.key == "*STEP") {
        model.step = Step{};
        if (kw.params.count("NAME")) model.step.name = kw.params.at("NAME");
        if (kw.params.count("NLGEOM") && upper(kw.params.at("NLGEOM")) == "YES") {
          model.step.type = AnalysisType::NonlinearStatic;
        }
        if (kw.params.count("AMPLITUDE")) model.step.amplitudeName = kw.params.at("AMPLITUDE");
      } else if (kw.key == "*PLASTIC") {
        plasticMode = kw.params.count("HARDENING") ? upper(kw.params.at("HARDENING")) : "BILINEAR";
      } else if (kw.key == "*DENSITY" || kw.key == "*EXPANSION" || kw.key == "*CONDUCTIVITY" || kw.key == "*SPECIFIC HEAT") {
        // 材料扩展参数在数据行读取。
      } else if (kw.key == "*AMPLITUDE") {
        currentAmplitude = kw.params.count("NAME") ? kw.params.at("NAME") : "AMP-1";
        model.amplitudes[currentAmplitude].name = currentAmplitude;
      } else if (kw.key == "*CONTACT PAIR") {
        ContactPair c;
        c.masterSurface = kw.params.count("MASTER") ? kw.params.at("MASTER") : "";
        c.slaveSurface = kw.params.count("SLAVE") ? kw.params.at("SLAVE") : "";
        if (kw.params.count("PENALTY")) c.penalty = std::stod(kw.params.at("PENALTY"));
        if (kw.params.count("FRICTION")) c.friction = std::stod(kw.params.at("FRICTION"));
        model.contacts.push_back(c);
      } else if (kw.key == "*CLOAD" || kw.key == "*DLOAD" || kw.key == "*DSLOAD") {
        currentLoadAmplitude = kw.params.count("AMPLITUDE") ? kw.params.at("AMPLITUDE") : "";
      } else if (kw.key == "*COUPLING") {
        Coupling c;
        if (kw.params.count("REF NODE")) c.referenceNode = std::stoi(kw.params.at("REF NODE"));
        if (kw.params.count("SURFACE")) c.surfaceName = kw.params.at("SURFACE");
        if (kw.params.count("PENALTY")) c.penalty = std::stod(kw.params.at("PENALTY"));
        model.couplings.push_back(c);
      } else if (kw.key == "*MPC") {
        // MPC 关键字本身仅切换解析状态，具体方程在随后的数据行读取。
      } else if (kw.key == "*NODE OUTPUT" || kw.key == "*ELEMENT OUTPUT") {
        OutputRequest req;
        req.kind = kw.key;
        model.outputRequests.push_back(req);
      } else if (kw.key == "*STATIC") {
        if (model.step.type != AnalysisType::NonlinearStatic) model.step.type = AnalysisType::LinearStatic;
        if (kw.params.count("RIKS")) {
          model.step.useArcLength = true;
          model.step.type = AnalysisType::NonlinearStatic;
        }
      }
      continue;
    }

    if (currentSection == "*HEADING") {
      if (model.modelName == "Model-1") model.modelName = line;
    } else if (currentSection == "*NODE") {
      const auto p = splitCSV(line);
      if (p.size() < 4) throw std::runtime_error("NODE line requires id,x,y,z");
      model.nodes.push_back(Node{std::stoi(p[0]), {std::stod(p[1]), std::stod(p[2]), std::stod(p[3])}});
    } else if (currentSection == "*NSET") {
      appendSetItems(model.nsets, currentNset, splitCSV(line));
    } else if (currentSection == "*ELSET") {
      appendSetItems(model.elsets, currentElset, splitCSV(line));
    } else if (currentSection == "*AMPLITUDE") {
      auto p = splitCSV(line);
      for (size_t i = 0; i + 1 < p.size(); i += 2) {
        model.amplitudes[currentAmplitude].points.push_back({std::stod(p[i]), std::stod(p[i + 1])});
      }
    } else if (currentSection == "*ELEMENT") {
      const auto p = splitCSV(line);
      if (p.size() < 3) throw std::runtime_error("ELEMENT line has insufficient node ids");
      Element e;
      e.id = std::stoi(p[0]);
      e.type = currentElementType;
      e.material = currentMaterial.empty() ? "MAT1" : currentMaterial;
      for (size_t i = 1; i < p.size(); ++i) e.conn.push_back(std::stoi(p[i]));
      if (e.type == ElementType::Solid8) {
        e.state.eqPlasticStrain.assign(1, 0.0);
        e.state.stress.assign(1, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
      }
      model.elements.push_back(e);
      model.elsets[currentElset].push_back(e.id);
    } else if (currentSection == "*TRUSS SECTION" || currentSection == "*BEAM SECTION") {
      auto p = splitCSV(line);
      if (!p.empty()) {
        const double area = std::stod(p[0]);
        for (auto& e : model.elements) {
          if (e.area == 1.0) e.area = area;
          if (currentSection == "*BEAM SECTION") e.type = ElementType::Beam2;
        }
      }
    } else if (currentSection == "*SHELL SECTION") {
      auto p = splitCSV(line);
      if (!p.empty()) {
        const double t = std::stod(p[0]);
        for (auto& e : model.elements)
          if (e.type == ElementType::Shell4) e.thickness = t;
      }
    } else if (currentSection == "*ELASTIC") {
      auto p = splitCSV(line);
      if (p.empty()) throw std::runtime_error("ELASTIC requires E");
      if (currentMaterial.empty()) currentMaterial = "MAT1";
      auto& m = model.materials[currentMaterial];
      m.name = currentMaterial;
      m.young = std::stod(p[0]);
      if (p.size() >= 2) m.poisson = std::stod(p[1]);
      m.law = MaterialLaw::LinearElastic;
    } else if (currentSection == "*PLASTIC") {
      auto p = splitCSV(line);
      auto& m = model.materials[currentMaterial];
      if (p.size() >= 2) {
        m.yieldStress = std::stod(p[0]);
        m.hardening = std::stod(p[1]);
      }
      if (p.size() >= 1) {
        const double stress = std::stod(p[0]);
        const double peeq = (p.size() >= 2 ? std::stod(p[1]) : 0.0);
        m.plasticTable.push_back({peeq, stress});
      }
      m.law = (plasticMode == "J2") ? MaterialLaw::J2Plasticity : MaterialLaw::BilinearElastoPlastic;
    } else if (currentSection == "*DENSITY") {
      auto p = splitCSV(line);
      if (!p.empty()) model.materials[currentMaterial].density = std::stod(p[0]);
    } else if (currentSection == "*EXPANSION") {
      auto p = splitCSV(line);
      if (!p.empty()) model.materials[currentMaterial].expansion = std::stod(p[0]);
    } else if (currentSection == "*CONDUCTIVITY") {
      auto p = splitCSV(line);
      if (!p.empty()) model.materials[currentMaterial].conductivity = std::stod(p[0]);
    } else if (currentSection == "*SPECIFIC HEAT") {
      auto p = splitCSV(line);
      if (!p.empty()) model.materials[currentMaterial].specificHeat = std::stod(p[0]);
    } else if (currentSection == "*TEMPERATURE") {
      auto p = splitCSV(line);
      if (p.size() >= 2) {
        if (isIntegerToken(p[0])) model.temperatureBCs.push_back({std::stoi(p[0]), "", std::stod(p[1])});
        else model.temperatureBCs.push_back({0, p[0], std::stod(p[1])});
      }
    } else if (currentSection == "*DFLUX") {
      auto p = splitCSV(line);
      if (p.size() >= 2) {
        if (isIntegerToken(p[0])) model.heatLoads.push_back({std::stoi(p[0]), "", std::stod(p.back())});
        else model.heatLoads.push_back({0, p[0], std::stod(p.back())});
      }
    } else if (currentSection == "*BOUNDARY") {
      auto p = splitCSV(line);
      if (p.size() < 3) throw std::runtime_error("BOUNDARY requires target,dofStart,dofEnd[,value]");
      int d1 = std::stoi(p[1]);
      int d2 = std::stoi(p[2]);
      double v = p.size() >= 4 ? std::stod(p[3]) : 0.0;
      for (int d = d1; d <= d2; ++d) {
        if (d < 1 || d > kDofPerNode) continue;
        if (isIntegerToken(p[0])) {
          model.bcs.push_back({std::stoi(p[0]), "", d, v});
        } else {
          model.bcs.push_back({0, p[0], d, v});
        }
      }
    } else if (currentSection == "*CLOAD") {
      auto p = splitCSV(line);
      if (p.size() < 3) throw std::runtime_error("CLOAD requires target,dof,value");
      int dof = std::stoi(p[1]);
      if (dof < 1 || dof > kDofPerNode) throw std::runtime_error("Unsupported dof in CLOAD");
      if (isIntegerToken(p[0])) {
        model.loads.push_back({std::stoi(p[0]), "", dof, std::stod(p[2]), currentLoadAmplitude});
      } else {
        model.loads.push_back({0, p[0], dof, std::stod(p[2]), currentLoadAmplitude});
      }

    } else if (currentSection == "*CONTACT PAIR" && !model.contacts.empty()) {
      auto p = splitCSV(line);
      if (p.size() >= 2) {
        if (model.contacts.back().slaveSurface.empty()) model.contacts.back().slaveSurface = p[0];
        if (model.contacts.back().masterSurface.empty()) model.contacts.back().masterSurface = p[1];
      }
    } else if (currentSection == "*DLOAD" || currentSection == "*DSLOAD") {
      auto p = splitCSV(line);
      if (p.size() >= 3) {
        BodyLoad b;
        b.elset = p[0];
        b.type = upper(p[1]);
        b.value = std::stod(p[2]);
        b.amplitudeName = currentLoadAmplitude;
        model.bodyLoads.push_back(b);
      }
    } else if (currentSection == "*KINEMATIC" && !model.couplings.empty()) {
      auto p = splitCSV(line);
      if (p.size() >= 2) {
        int d1 = std::stoi(p[0]), d2 = std::stoi(p[1]);
        for (int d = d1; d <= d2; ++d) model.couplings.back().dofs.push_back(d);
      }
    } else if (currentSection == "*MPC") {
      // 读取一条 MPC 关系：slave = coef*master + offset（当前实现单主点扩展格式）。
      auto p = splitCSV(line);
      if (p.size() >= 4) {
        MpcConstraint m;
        m.slaveNode = std::stoi(p[0]);
        m.slaveDof = std::stoi(p[1]);
        m.masterNodes.push_back(std::stoi(p[2]));
        m.masterDofs.push_back(std::stoi(p[3]));
        m.coefficients.push_back(p.size() >= 5 ? std::stod(p[4]) : 1.0);
        if (p.size() >= 6) m.offset = std::stod(p[5]);
        if (m.slaveDof >= 1 && m.slaveDof <= kDofPerNode && m.masterDofs[0] >= 1 && m.masterDofs[0] <= kDofPerNode) {
          model.mpcs.push_back(m);
        }
      }
    } else if ((currentSection == "*NODE OUTPUT" || currentSection == "*ELEMENT OUTPUT") &&
               !model.outputRequests.empty()) {
      auto p = splitCSV(line);
      for (const auto& v : p)
        if (!v.empty()) model.outputRequests.back().variables.push_back(upper(v));
    } else if (currentSection == "*STATIC") {
      auto p = splitCSV(line);
      if (!p.empty()) {
        model.step.initialIncrement = std::max(1e-12, std::stod(p[0]));
        model.step.increments = std::max(1, static_cast<int>(1.0 / model.step.initialIncrement));
      }
      if (p.size() >= 2) model.step.totalLoadFactor = std::max(1e-12, std::stod(p[1]));
      if (p.size() >= 3) model.step.minIncrement = std::max(1e-12, std::stod(p[2]));
      if (p.size() >= 4) model.step.maxIncrement = std::max(model.step.minIncrement, std::stod(p[3]));
      if (model.step.useArcLength && p.size() >= 2) {
        model.step.arcLengthRadius = std::max(1e-12, std::stod(p[1]));
      }
    } else if (currentSection == "*CONTROLS") {
      // 扩展控制参数：Newton 迭代、容差、cutback、弧长半径调节参数。
      auto p = splitCSV(line);
      if (p.size() >= 2) {
        model.step.maxNewtonIters = std::stoi(p[0]);
        model.step.tolerance = std::stod(p[1]);
      }
      if (p.size() >= 3) model.step.maxCutbacks = std::stoi(p[2]);
      if (p.size() >= 4) model.step.arcLengthGrowFactor = std::stod(p[3]);
      if (p.size() >= 5) model.step.arcLengthShrinkFactor = std::stod(p[4]);
      if (p.size() >= 6) model.step.arcLengthMinRadius = std::stod(p[5]);
      if (p.size() >= 7) model.step.arcLengthMaxRadius = std::stod(p[6]);
    }
  }

  if (model.materials.empty()) model.materials["MAT1"] = Material{"MAT1"};
  for (auto& c : model.couplings) {
    if (!c.surfaceName.empty()) {
      auto it = model.nsets.find(c.surfaceName);
      if (it != model.nsets.end()) c.surfaceNodes = it->second;
    }
  }
  for (auto& e : model.elements) {
    if (e.material.empty()) e.material = model.materials.begin()->first;
    for (int nid : e.conn) {
      if (nodeIndex(model, nid) < 0) throw std::runtime_error("Element references missing node id");
    }
  }
  for (const auto& m : model.mpcs) {
    if (nodeIndex(model, m.slaveNode) < 0) throw std::runtime_error("MPC references missing slave node id");
    for (int mn : m.masterNodes) {
      if (nodeIndex(model, mn) < 0) throw std::runtime_error("MPC references missing master node id");
    }
  }

  return model;
}

}  // namespace fem
