from amalgkit.arg_utils import clone_namespace, namespace_to_dict


class TestArgUtils:
    def test_namespace_to_dict_includes_class_attributes(self):
        class Args:
            out_dir = '/tmp/out'
            threads = 4

        data = namespace_to_dict(Args())

        assert data['out_dir'] == '/tmp/out'
        assert data['threads'] == 4

    def test_clone_namespace_overrides_copied_attributes(self):
        class Args:
            out_dir = '/tmp/out'
            threads = 4

        cloned = clone_namespace(Args(), threads=1, internal_jobs=2)

        assert cloned.out_dir == '/tmp/out'
        assert cloned.threads == 1
        assert cloned.internal_jobs == 2
